//! ASSIST propagator adapter.
//!
//! Implements the [`Propagator`] trait using the `assist-rs` crate
//! (ASSIST + REBOUND N-body integration). First-order machinery only:
//! STMs from `assist_propagate` drive covariance propagation via
//! `Σ(t) = Φ·Σ₀·Φᵀ` (so every covariance is tagged
//! [`CovarianceKind::Linear`]), and the observation Jacobian is composed
//! as `J_local · R_ecl→eq · Φ` when requested. Sample-based covariance
//! methods (SigmaPoint / MonteCarlo) are loud
//! [`PropagatorError::Unsupported`] errors, as is the Fast integrator
//! profile — IAS15 is the only integrator.
//!
//! Gated behind the `assist` feature.

use std::path::Path;

use assist_rs::coordinates::cartesian_to_spherical_jacobian;
use assist_rs::{
    AssistData, Ephemeris as AssistEphemeris, Ias15AdaptiveMode, IntegratorConfig,
    ObservatoryTable, Observer as AssistObserver, Orbit as AssistOrbit, Origin as AssistOrigin,
    assist_generate_ephemeris_single, assist_get_state, assist_propagate_single,
    ecliptic_to_equatorial, equatorial_to_ecliptic,
};

use crate::linalg::RAD2DEG;
use crate::linalg::{mat6_min_eigenvalue_sym, propagate_covariance_6x6};
use thor_rs_propagator::PropagatorError;
use thor_rs_propagator::{
    CartesianState, CovarianceKind, CovarianceMethod, CovarianceQuality, Ephemeris,
    EphemerisConfig, Frame, IntegratorProfile, ObserverState, Origin, PropagatedState,
    PropagationConfig, PropagationProfile, Propagator, SphericalState, TestOrbit,
};

/// ASSIST propagator: N-body integration via assist-rs + REBOUND.
pub struct AssistPropagator {
    data: AssistData,
    integrator: IntegratorConfig,
}

/// Bound IAS15 timestep + precision so a single pathological LM update can't
/// drag the integrator into endless predictor-corrector retries. Empirically,
/// THOR's Phase-2 DC blew Phase 2 walltime from ~ms to minutes per orbit
/// because the unbounded adaptive-step path could shrink `dt` toward zero on
/// stiff candidate states. `min_dt = 1e-6 d` (~0.1 s) is well below the IAS15
/// step a healthy 30-day arc takes, but bounds the worst case to a finite
/// number of steps. `epsilon = 1e-9` matches REBOUND's default. `Prs23` is
/// REBOUND's recommended adaptive criterion since 2024-01.
///
/// Tried matching `adam-assist`'s looser pairing (`epsilon = 1e-6`,
/// `min_dt = 1e-9`, `adaptive_mode = Global`) to sidestep the silent
/// floor-hit case, but the larger steps push debug-mode propagation past
/// the assist ephemeris coverage and SIGSEGV `test_propagator_consistency`.
/// Tight defaults retained until validation can confirm parity.
fn default_integrator_config() -> IntegratorConfig {
    IntegratorConfig {
        initial_dt: Some(1e-6),
        min_dt: Some(1e-6),
        epsilon: Some(1e-9),
        adaptive_mode: Some(Ias15AdaptiveMode::Prs23),
    }
}

impl AssistPropagator {
    /// Build from a pre-loaded ASSIST ephemeris.
    ///
    /// The observatory table is unset; calls requiring observatory codes
    /// will fail until one is provided via [`AssistPropagator::with_obs_table`].
    pub fn new(ephem: AssistEphemeris) -> Self {
        Self {
            data: AssistData::new(ephem),
            integrator: default_integrator_config(),
        }
    }

    /// Attach an observatory table (required for observer lookups by code).
    pub fn with_obs_table(mut self, table: ObservatoryTable) -> Self {
        self.data = self.data.with_observatory(table);
        self
    }

    /// Override the IAS15 integrator configuration. Use sparingly — the
    /// default (`default_integrator_config`) is tuned to bound runaway
    /// adaptive-step shrinkage on stiff orbits.
    pub fn with_integrator_config(mut self, integrator: IntegratorConfig) -> Self {
        self.integrator = integrator;
        self
    }

    /// Load from two SPK files (`de440.bsp` + `sb441-n16.bsp`).
    pub fn from_paths(
        planets_path: impl AsRef<Path>,
        asteroids_path: impl AsRef<Path>,
    ) -> Result<Self, PropagatorError> {
        let ephem = AssistEphemeris::from_paths(planets_path.as_ref(), asteroids_path.as_ref())
            .map_err(|e| PropagatorError::Other(format!("Failed to load ASSIST ephemeris: {e}")))?;
        Ok(Self::new(ephem))
    }

    /// Load ephemeris and observatory table together. `eop_paths` is the
    /// three Earth-orientation PCK kernels (predict, historical, current
    /// high-precision) — without them, observer lookups for ground sites
    /// like LSST X05 fail per assist-rs's strict EOP requirement.
    pub fn from_paths_with_obs<P: AsRef<Path>>(
        planets_path: impl AsRef<Path>,
        asteroids_path: impl AsRef<Path>,
        obs_codes_path: impl AsRef<Path>,
        eop_paths: &[P],
    ) -> Result<Self, PropagatorError> {
        let prop = Self::from_paths(planets_path, asteroids_path)?;
        let table = ObservatoryTable::from_json(obs_codes_path.as_ref())
            .map_err(|e| PropagatorError::Other(format!("Failed to load obs codes: {e}")))?;
        let eo = assist_rs::earth_orientation::EarthOrientation::from_paths(eop_paths)
            .map_err(|e| PropagatorError::Other(format!("Failed to load EOP kernels: {e}")))?;
        Ok(prop.with_obs_table(table.with_earth_orientation(eo)))
    }

    /// Load everything via [`assist_rs::data::DataManager`], downloading any
    /// missing files to its default cache directory.
    pub fn from_data_manager() -> Result<Self, PropagatorError> {
        let dm = assist_rs::data::DataManager::new();
        let paths = dm
            .ensure_ready()
            .map_err(|e| PropagatorError::Other(format!("DataManager: {e}")))?;
        let eop = paths.eop_kernels();
        Self::from_paths_with_obs(&paths.planets, &paths.asteroids, &paths.obscodes, &eop)
    }

    /// Borrow the underlying ASSIST ephemeris.
    pub fn ephem(&self) -> &AssistEphemeris {
        &self.data.ephem
    }

    fn to_assist_orbit(orbit: &TestOrbit) -> Result<AssistOrbit, PropagatorError> {
        let s = &orbit.state;
        if s.frame != Frame::EclipticJ2000 {
            return Err(PropagatorError::InvalidOrbit(format!(
                "AssistPropagator requires EclipticJ2000 frame, got {:?}",
                s.frame
            )));
        }
        if s.origin != Origin::Sun {
            return Err(PropagatorError::InvalidOrbit(format!(
                "AssistPropagator requires Sun origin (heliocentric), got {:?}",
                s.origin
            )));
        }
        Ok(AssistOrbit::new([s.x, s.y, s.z, s.vx, s.vy, s.vz], s.epoch))
    }

    fn from_helio_ecliptic(state: [f64; 6], epoch: f64) -> CartesianState {
        CartesianState {
            x: state[0],
            y: state[1],
            z: state[2],
            vx: state[3],
            vy: state[4],
            vz: state[5],
            epoch,
            frame: Frame::EclipticJ2000,
            origin: Origin::Sun,
            covariance: None,
        }
    }

    /// Heliocentric ecliptic J2000 state of an origin center at `epoch`.
    /// Sun is identically zero; SSB and Earth are looked up from the
    /// ephemeris. Used to compose origin translations in
    /// [`Propagator::transform_state`] with the Sun as hub.
    fn origin_helio_ecl_state(
        &self,
        origin: Origin,
        epoch: f64,
    ) -> Result<[f64; 6], PropagatorError> {
        let a_origin = match origin {
            Origin::Sun => return Ok([0.0; 6]),
            Origin::SolarSystemBarycenter => AssistOrigin::SolarSystemBarycenter,
            Origin::Earth => AssistOrigin::Earth,
        };
        let states = assist_get_state(&self.data, &a_origin, &[epoch], Some(1))
            .map_err(|e| PropagatorError::Other(format!("assist_get_state({origin:?}): {e}")))?;
        states
            .first()
            .map(|s| s.state)
            .ok_or_else(|| PropagatorError::Other(format!("no state returned for {origin:?}")))
    }
}

/// Validate a [`PropagationProfile`] against ASSIST's capabilities.
///
/// - `events`: documented no-op — ASSIST runs no event detection, so
///   both profile settings already describe its behavior.
/// - `integrator: Fast`: ASSIST has only IAS15. This errors loudly
///   unless the caller explicitly set `allow_accurate_substitute`, in
///   which case the substitution is logged once per process so the run
///   record reflects what actually integrated.
/// - `cache_dense`: documented no-op — assist-rs exposes no dense-step
///   cache.
fn check_profile(profile: &PropagationProfile) -> Result<(), PropagatorError> {
    if profile.integrator == IntegratorProfile::Fast {
        if !profile.allow_accurate_substitute {
            return Err(PropagatorError::Unsupported(
                "AssistPropagator has no Fast (survey-grade) integrator — only IAS15. \
                 Set PropagationProfile::allow_accurate_substitute to run IAS15 instead, \
                 or use the empyrean backend for DOP853."
                    .to_string(),
            ));
        }
        static SUBSTITUTE_LOGGED: std::sync::atomic::AtomicBool =
            std::sync::atomic::AtomicBool::new(false);
        if !SUBSTITUTE_LOGGED.swap(true, std::sync::atomic::Ordering::Relaxed) {
            log::warn!(
                "PropagationProfile requested the Fast integrator; ASSIST substitutes IAS15 \
                 (allow_accurate_substitute set). Reported once per process."
            );
        }
    }
    Ok(())
}

/// Definiteness tag for a propagated covariance. Φ Σ₀ Φᵀ with PSD Σ₀ is
/// PSD up to round-off, so anything beyond a relative round-off tolerance
/// is tagged honestly as indefinite rather than silently passed on.
fn covariance_quality_6x6(cov: &[[f64; 6]; 6]) -> CovarianceQuality {
    let min_eig = mat6_min_eigenvalue_sym(cov);
    let max_diag = (0..6).map(|i| cov[i][i]).fold(0.0_f64, f64::max);
    if min_eig >= -1e-12 * max_diag.max(f64::MIN_POSITIVE) {
        CovarianceQuality::PositiveDefinite
    } else {
        CovarianceQuality::Indefinite { min_eig }
    }
}

impl Propagator for AssistPropagator {
    fn force_model(&self) -> &'static str {
        // Reproducibility provenance: identifies the integrator backend
        // and the deterministic configuration used. Update when the
        // ephemeris or integrator scheme changes.
        "assist:de440+sb441-n16+ias15-prs23"
    }

    fn propagate(
        &self,
        orbit: &TestOrbit,
        epochs: &[f64],
        config: &PropagationConfig,
    ) -> Result<Vec<PropagatedState>, PropagatorError> {
        check_profile(&config.profile)?;
        let a_orbit = Self::to_assist_orbit(orbit)?;

        // Request STMs whenever we need them: either the caller asked or the
        // input has covariance that we need to propagate forward.
        let want_stm = config.compute_stm || orbit.state.covariance.is_some();
        let results =
            assist_propagate_single(&self.data, &a_orbit, epochs, want_stm, &self.integrator)
                .map_err(|e| PropagatorError::PropagationFailed(format!("{e}")))?;

        let mut out = Vec::with_capacity(results.len());
        for r in results.into_iter() {
            let cov = match (orbit.state.covariance, r.stm.as_ref()) {
                (Some(c0), Some(phi)) => Some(propagate_covariance_6x6(phi, &c0)),
                _ => None,
            };
            let stm = if config.compute_stm { r.stm } else { None };
            let mut state = Self::from_helio_ecliptic(r.state, r.epoch);
            state.covariance = cov;
            // Honest provenance tags: ASSIST's only machinery is the linear
            // Φ Σ₀ Φᵀ mapping, verified for definiteness before tagging.
            let (covariance_kind, covariance_quality) = match &cov {
                Some(c) => (
                    Some(CovarianceKind::Linear),
                    Some(covariance_quality_6x6(c)),
                ),
                None => (None, None),
            };
            out.push(PropagatedState {
                state,
                stm,
                covariance_kind,
                covariance_quality,
            });
        }
        Ok(out)
    }

    fn compute_observers(
        &self,
        codes: &[String],
        epochs: &[f64],
    ) -> Result<Vec<ObserverState>, PropagatorError> {
        if codes.len() != epochs.len() {
            return Err(PropagatorError::Other(format!(
                "compute_observers: codes.len()={} != epochs.len()={}",
                codes.len(),
                epochs.len()
            )));
        }

        let mut out = Vec::with_capacity(codes.len());
        for (code, &epoch) in codes.iter().zip(epochs.iter()) {
            let origin = AssistOrigin::Observatory(code.clone());
            let states = assist_get_state(&self.data, &origin, &[epoch], Some(1))
                .map_err(|e| PropagatorError::InvalidObserver(format!("{code}: {e}")))?;
            let s = states
                .into_iter()
                .next()
                .ok_or_else(|| PropagatorError::InvalidObserver(format!("no state for {code}")))?;
            out.push(ObserverState {
                code: code.clone(),
                state: Self::from_helio_ecliptic(s.state, s.epoch),
            });
        }
        Ok(out)
    }

    fn compute_body_positions(
        &self,
        body: &str,
        epochs: &[f64],
    ) -> Result<Vec<[f64; 3]>, PropagatorError> {
        let origin = match body.to_lowercase().as_str() {
            "sun" => AssistOrigin::Sun,
            "mercury" => AssistOrigin::MercuryBarycenter,
            "venus" => AssistOrigin::VenusBarycenter,
            "earth" => AssistOrigin::Earth,
            "moon" => AssistOrigin::Moon,
            "mars" => AssistOrigin::MarsBarycenter,
            "jupiter" => AssistOrigin::JupiterBarycenter,
            "saturn" => AssistOrigin::SaturnBarycenter,
            "uranus" => AssistOrigin::UranusBarycenter,
            "neptune" => AssistOrigin::NeptuneBarycenter,
            "pluto" => AssistOrigin::PlutoBarycenter,
            other => {
                return Err(PropagatorError::Other(format!(
                    "compute_body_positions: unknown body '{other}' (expected sun/mercury/venus/earth/moon/mars/jupiter/saturn/uranus/neptune/pluto)"
                )));
            }
        };
        let states = assist_get_state(&self.data, &origin, epochs, Some(epochs.len()))
            .map_err(|e| PropagatorError::Other(format!("compute_body_positions({body}): {e}")))?;
        if states.len() != epochs.len() {
            return Err(PropagatorError::Other(format!(
                "compute_body_positions({body}): expected {} states, got {}",
                epochs.len(),
                states.len()
            )));
        }
        Ok(states
            .into_iter()
            .map(|s| [s.state[0], s.state[1], s.state[2]])
            .collect())
    }

    fn transform_state(
        &self,
        state: &CartesianState,
        target_frame: Frame,
        target_origin: Origin,
    ) -> Result<CartesianState, PropagatorError> {
        let mut x = [state.x, state.y, state.z, state.vx, state.vy, state.vz];
        // Covariance carries through the whole transform: identity and pure
        // origin translations leave it untouched; frame rotations apply the
        // exact congruence R Σ Rᵀ. (Previously every path — including the
        // identity no-op — silently stripped it, which is exactly the path
        // covariant states traverse via on_sky_cov's normalisation calls.)
        let mut cov = state.covariance;
        let mut cur_frame = state.frame;
        let mut cur_origin = state.origin;

        // Rotate frame first: pure geometric rotation, does not touch origin.
        // Exhaustive over (from, to) so a new Frame variant is a compile
        // error here instead of a latent runtime panic; the same-frame arms
        // are dead under the `!=` guard.
        if cur_frame != target_frame {
            match (cur_frame, target_frame) {
                (Frame::EclipticJ2000, Frame::Equatorial) => {
                    x = ecliptic_to_equatorial(&x);
                    cov = cov.map(|c| propagate_covariance_6x6(&ecl_to_eq_6x6(), &c));
                }
                (Frame::Equatorial, Frame::EclipticJ2000) => {
                    x = equatorial_to_ecliptic(&x);
                    // The inverse of an orthonormal rotation is its transpose.
                    cov = cov.map(|c| {
                        propagate_covariance_6x6(
                            &hyperjet::linalg::mat6::mat6_transpose(&ecl_to_eq_6x6()),
                            &c,
                        )
                    });
                }
                (Frame::EclipticJ2000, Frame::EclipticJ2000)
                | (Frame::Equatorial, Frame::Equatorial) => {}
            }
            cur_frame = target_frame;
        }

        // Shift origin: compose through the Sun as hub (Sun's own state is
        // identically zero, SSB and Earth come from the ephemeris), so every
        // (from, to) pair among Sun/SSB/Earth is handled. A pure translation
        // leaves the covariance untouched.
        if cur_origin != target_origin {
            let cur_helio = self.origin_helio_ecl_state(cur_origin, state.epoch)?;
            let tgt_helio = self.origin_helio_ecl_state(target_origin, state.epoch)?;
            // Offset = new_origin_pos − old_origin_pos, in heliocentric ecliptic.
            let mut offset_ecl = [0.0_f64; 6];
            for i in 0..6 {
                offset_ecl[i] = tgt_helio[i] - cur_helio[i];
            }
            let offset = if cur_frame == Frame::EclipticJ2000 {
                offset_ecl
            } else {
                ecliptic_to_equatorial(&offset_ecl)
            };
            for i in 0..6 {
                x[i] -= offset[i];
            }
            cur_origin = target_origin;
        }

        Ok(CartesianState {
            x: x[0],
            y: x[1],
            z: x[2],
            vx: x[3],
            vy: x[4],
            vz: x[5],
            epoch: state.epoch,
            frame: cur_frame,
            origin: cur_origin,
            covariance: cov,
        })
    }

    fn generate_ephemeris(
        &self,
        orbit: &TestOrbit,
        observers: &[ObserverState],
        config: &EphemerisConfig,
    ) -> Result<Vec<Ephemeris>, PropagatorError> {
        // Covariance output is a real config contract (never keyed off the
        // orbit payload alone): `compute_covariance = true` requires an input
        // covariance and a method ASSIST actually has; `false` skips the
        // covariance pass entirely, even for covariant orbits.
        let want_cov = config.compute_covariance;
        let want_jac = config.compute_jacobian;
        if want_cov {
            if orbit.state.covariance.is_none() {
                return Err(PropagatorError::MissingCovariance(format!(
                    "AssistPropagator::generate_ephemeris: compute_covariance requested \
                     for orbit {} but the orbit carries no covariance",
                    orbit.id
                )));
            }
            match config.covariance_method {
                // Analytic is ASSIST's native (and only) method; Auto resolves
                // to the same first-order machinery — an explicit capability
                // statement (see CovarianceMethod docs), not a silent downgrade.
                CovarianceMethod::Auto | CovarianceMethod::Analytic => {}
                CovarianceMethod::SigmaPoint | CovarianceMethod::MonteCarlo => {
                    return Err(PropagatorError::Unsupported(format!(
                        "AssistPropagator has no {:?} covariance machinery — only the \
                         first-order STM (Analytic). Use the empyrean backend for \
                         sample-based covariance propagation.",
                        config.covariance_method
                    )));
                }
            }
        }

        let a_orbit = Self::to_assist_orbit(orbit)?;

        let a_observers: Vec<AssistObserver> = observers
            .iter()
            .map(|obs| {
                AssistObserver::new(AssistOrigin::Observatory(obs.code.clone()), obs.state.epoch)
            })
            .collect();

        let results = assist_generate_ephemeris_single(
            &self.data,
            &a_orbit,
            &a_observers,
            Some(1),
            &self.integrator,
        )
        .map_err(|e| PropagatorError::PropagationFailed(format!("{e}")))?;

        // STM at t_emit is needed for both covariance propagation
        // (Σ(t_emit) = Φ(t_emit, t₀)·Σ₀·Φᵀ) and observation Jacobian
        // (∂spherical/∂x₀ = J_local · R_ecl→eq · Φ). t_emit = t_obs − τ.
        // Using t_obs would attach Σ(t_obs) to a state at t_emit, introducing
        // an O(τ·dynamics) discrepancy — negligible for MBA (τ ≈ 0.01 d) but
        // ~1e-3 fractional for TNO (τ ≈ 0.2 d). One propagate call serves
        // both consumers below.
        let stms: Vec<Option<[[f64; 6]; 6]>> = if want_cov || want_jac {
            let prop_epochs: Vec<f64> = observers
                .iter()
                .zip(results.iter())
                .map(|(o, r)| o.state.epoch - r.light_time)
                .collect();
            let prop =
                assist_propagate_single(&self.data, &a_orbit, &prop_epochs, true, &self.integrator)
                    .map_err(|e| PropagatorError::PropagationFailed(format!("{e}")))?;
            prop.into_iter().map(|p| p.stm).collect()
        } else {
            vec![None; results.len()]
        };
        // Under the covariance contract a missing STM cannot silently become
        // a missing covariance — that is a backend failure, reported loudly.
        if want_cov && stms.iter().any(|s| s.is_none()) {
            return Err(PropagatorError::PropagationFailed(
                "AssistPropagator::generate_ephemeris: compute_covariance requested but \
                 assist_propagate returned no STM for one or more emission epochs"
                    .to_string(),
            ));
        }

        let prop_covs: Vec<Option<[[f64; 6]; 6]>> = if want_cov {
            let cov0 = orbit
                .state
                .covariance
                .expect("checked Some above under want_cov");
            stms.iter()
                .map(|phi| phi.map(|p| propagate_covariance_6x6(&p, &cov0)))
                .collect()
        } else {
            vec![None; results.len()]
        };

        // The observation Jacobian blocks serve two consumers: the DC
        // (`compute_jacobian`) and the on-sky spherical covariance
        // Σ_sph = J_obs Σ₀ J_obsᵀ (`compute_covariance`) — the marginal cost
        // over the STM already in hand is one 6×6 triple product per epoch.
        let obs_jacs: Vec<Option<[[f64; 6]; 6]>> = if want_jac || want_cov {
            stms.iter()
                .zip(results.iter())
                .zip(observers.iter())
                .map(|((phi_opt, r), obs)| {
                    phi_opt.map(|phi| {
                        observation_jacobian_assist(&phi, &r.aberrated_state, &obs.state)
                    })
                })
                .collect()
        } else {
            vec![None; results.len()]
        };

        let mut out = Vec::with_capacity(results.len());
        for (i, (r, obs)) in results.iter().zip(observers.iter()).enumerate() {
            // assist-rs spherical order: [rho AU, ra rad, dec rad, drho, dra, ddec].
            // THOR SphericalState stores angles in degrees.
            // On-sky spherical covariance (R1): J_obs maps Σ₀ directly into
            // the spherical observable coordinates (deg² angular blocks, raw
            // RA — the SphericalState convention). This is what the
            // Mahalanobis observation filter and sky-ellipse pre-prune
            // consume; without it they silently degrade to a bare cone.
            let sph_cov = if want_cov {
                let j_obs = obs_jacs[i].expect("STM presence checked above under want_cov");
                let cov0 = orbit
                    .state
                    .covariance
                    .expect("checked Some above under want_cov");
                Some(propagate_covariance_6x6(&j_obs, &cov0))
            } else {
                None
            };
            let sph = SphericalState {
                rho: r.spherical[0],
                lon: r.spherical[1] * RAD2DEG,
                lat: r.spherical[2] * RAD2DEG,
                vrho: r.spherical[3],
                vlon: r.spherical[4] * RAD2DEG,
                vlat: r.spherical[5] * RAD2DEG,
                epoch: r.epoch,
                frame: Frame::Equatorial,
                origin: Origin::Sun,
                covariance: sph_cov,
            };
            // Aberrated state has its position/velocity at the emission
            // epoch (light-time corrected), so label it `t_emit` for
            // internal self-consistency. Consumers that re-propagate from
            // this state (existed in earlier branches; landmine for
            // future code) thus start at the right epoch.
            let t_emit = r.epoch - r.light_time;
            let mut aberrated = Self::from_helio_ecliptic(r.aberrated_state, t_emit);
            aberrated.covariance = prop_covs[i];
            out.push(Ephemeris {
                state: sph,
                aberrated_state: aberrated,
                observer_state: obs.state,
                light_time: Some(r.light_time),
                // Only surfaced when the DC asked for it; the covariance-only
                // path uses the block internally for Σ_sph and keeps the
                // output slot None so consumers see exactly what they requested.
                observation_jacobian: if want_jac { obs_jacs[i] } else { None },
                // ASSIST has no second-order (STT) machinery, so there is no
                // honest mean shift to report. None — which the usable-arc
                // nonlinearity gate treats as zero shift; arc.rs logs the
                // degraded gate once per run.
                mean_shift: None,
            });
        }
        Ok(out)
    }
}

/// Build the 6×6 observation Jacobian
/// `J_obs = ∂(ρ, α, δ, ρ̇, α̇, δ̇)_t_emit / ∂x_helio_ecl_t₀`.
///
/// Composition (left-to-right):
/// 1. `J_local`: cartesian-to-spherical Jacobian at the equatorial topocentric
///    state — from `assist-rs`. Output rows: spherical (rad). Input columns:
///    equatorial topocentric `(x, y, z, vx, vy, vz)`.
/// 2. `R_ecl→eq`: 6×6 block-diagonal rotation. ASSIST exposes the vector form
///    `ecliptic_to_equatorial`; the matrix is reconstructed here from the
///    canonical mean obliquity ε.
/// 3. `Φ`: state transition matrix at `t_emit`, heliocentric ecliptic.
///
/// Final rows for angular outputs (RA, Dec, RA-rate, Dec-rate) are scaled
/// `rad → deg` so the Jacobian's units match `Ephemeris::state.lon/lat`
/// (degrees) — consumers take `J[1][k]` directly as `∂α_deg/∂x₀[k]`.
fn observation_jacobian_assist(
    phi: &[[f64; 6]; 6],
    aberrated_helio_ecl: &[f64; 6],
    observer: &CartesianState,
) -> [[f64; 6]; 6] {
    let aberrated_eq = ecliptic_to_equatorial(aberrated_helio_ecl);
    let observer_ecl = [
        observer.x,
        observer.y,
        observer.z,
        observer.vx,
        observer.vy,
        observer.vz,
    ];
    let observer_eq = ecliptic_to_equatorial(&observer_ecl);

    let dx = [
        aberrated_eq[0] - observer_eq[0],
        aberrated_eq[1] - observer_eq[1],
        aberrated_eq[2] - observer_eq[2],
    ];
    let dv = [
        aberrated_eq[3] - observer_eq[3],
        aberrated_eq[4] - observer_eq[4],
        aberrated_eq[5] - observer_eq[5],
    ];
    let j_local = cartesian_to_spherical_jacobian(dx, dv);

    let r = ecl_to_eq_6x6();
    let tmp = mat6_mul(&j_local, &r);
    let mut j_obs = mat6_mul(&tmp, phi);

    // Convert angular rows (RA, Dec, RA-rate, Dec-rate) from rad to deg.
    for k in 0..6 {
        j_obs[1][k] *= RAD2DEG;
        j_obs[2][k] *= RAD2DEG;
        j_obs[4][k] *= RAD2DEG;
        j_obs[5][k] *= RAD2DEG;
    }
    j_obs
}

/// 6×6 block-diagonal rotation taking heliocentric ecliptic J2000 vectors to
/// equatorial. Inverse of `assist_rs::coordinates::equatorial_to_ecliptic`
/// in matrix form (rotates about x-axis by +ε).
fn ecl_to_eq_6x6() -> [[f64; 6]; 6] {
    // Use the same constants as assist-rs (J2000 mean obliquity ε ≈ 23.4393°).
    const COS_EPS: f64 = 0.917_482_062_070_108_2;
    const SIN_EPS: f64 = 0.397_777_155_929_776_9;
    let mut r = [[0.0f64; 6]; 6];
    r[0][0] = 1.0;
    r[1][1] = COS_EPS;
    r[1][2] = -SIN_EPS;
    r[2][1] = SIN_EPS;
    r[2][2] = COS_EPS;
    r[3][3] = 1.0;
    r[4][4] = COS_EPS;
    r[4][5] = -SIN_EPS;
    r[5][4] = SIN_EPS;
    r[5][5] = COS_EPS;
    r
}

#[inline]
fn mat6_mul(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for k in 0..6 {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
    c
}
