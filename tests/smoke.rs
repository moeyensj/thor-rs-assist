//! Construction + propagation smoke test against real ephemeris data.
//!
//! Needs the assist-rs data directory (~/.cache/assist-rs); CI downloads the
//! SPICE kernels and writes a minimal geocentric obscodes table. Run
//! single-threaded if extended: REBOUND's SPK reader is not thread-safe.

use thor_rs_propagator::{
    CartesianState, EphemerisConfig, Frame, Origin, PropagationConfig, Propagator, TestOrbit,
};

fn mba() -> TestOrbit {
    TestOrbit {
        id: "smoke_mba".into(),
        object_id: None,
        bundle_id: None,
        nside: 0,
        state: CartesianState {
            x: 2.3,
            y: 1.0,
            z: 0.1,
            vx: -0.005,
            vy: 0.009,
            vz: 0.0001,
            epoch: 60800.0,
            frame: Frame::EclipticJ2000,
            origin: Origin::Sun,
            covariance: None,
        },
    }
}

#[test]
fn propagate_and_ephemeris_roundtrip() {
    let prop = thor_rs_assist::AssistPropagator::from_data_manager()
        .expect("assist data (kernels + obscodes) must be present");
    let states = prop
        .propagate(&mba(), &[60810.0], &PropagationConfig::default())
        .expect("propagate");
    assert_eq!(states.len(), 1);
    let s = &states[0].state;
    assert!(s.x.is_finite() && s.y.is_finite() && s.z.is_finite());
    let r = (s.x * s.x + s.y * s.y + s.z * s.z).sqrt();
    assert!((1.5..4.0).contains(&r), "MBA stayed in the belt: r = {r}");

    let observers = prop
        .compute_observers(&["500".into()], &[60810.0])
        .expect("observers");
    let eph = prop
        .generate_ephemeris(&mba(), &observers, &EphemerisConfig::default())
        .expect("ephemeris");
    assert_eq!(eph.len(), 1);
    assert!(eph[0].state.lon.is_finite() && eph[0].state.lat.is_finite());
    assert_eq!(prop.force_model(), "assist:de440+sb441-n16+ias15-prs23");
}
