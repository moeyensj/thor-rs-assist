//! thor-rs-assist: the fully-open ASSIST backend for THOR.
//!
//! Implements [`thor_rs_propagator::Propagator`] over
//! [ASSIST](https://github.com/matthewholman/assist) +
//! [REBOUND](https://github.com/hannorein/rebound) via
//! [assist-rs](https://github.com/B612-Asteroid-Institute/assist-rs):
//! IAS15 N-body propagation with finite-difference state transition
//! matrices, light-time-corrected astrometric ephemerides with analytic
//! observation Jacobians, and MPC observatory handling through assist-rs's
//! data manager.
//!
//! Licensing: this crate's own source is BSD-3-Clause, but assist-rs and the
//! ASSIST/REBOUND stack beneath it are GPL-3.0 — any binary that includes
//! this backend is a GPL-3.0 combined work, which distributors must account
//! for. THOR's release configuration keeps it out of prebuilt artifacts.

// Inherited from the THOR source this adapter was extracted from: matrix
// math uses index-based loops for clarity, and the trait methods carry full
// parameter sets.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

mod linalg;
mod propagator;

pub use propagator::AssistPropagator;
