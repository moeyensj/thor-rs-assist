# thor-rs-assist

The fully-open ASSIST/REBOUND backend for [THOR](https://github.com/moeyensj/thor) — IAS15 N-body propagation and astrometric ephemerides behind the [thor-rs-propagator](https://github.com/moeyensj/thor-rs-propagator) contract

<a href="https://github.com/moeyensj/thor-rs-assist/actions/workflows/rust.yml"><img src="https://github.com/moeyensj/thor-rs-assist/actions/workflows/rust.yml/badge.svg" alt="CI"></a>
<a href="https://crates.io/crates/thor-rs-assist"><img src="https://img.shields.io/crates/v/thor-rs-assist.svg?style=flat-square&label=crates.io" alt="crates.io"></a>
<a href="https://docs.rs/thor-rs-assist"><img src="https://img.shields.io/docsrs/thor-rs-assist?style=flat-square&label=docs.rs" alt="docs.rs"></a>
<br>
<a href="Cargo.toml"><img src="https://img.shields.io/badge/rustc-1.94%2B-orange?style=flat-square&logo=rust" alt="MSRV 1.94"></a>
<a href="LICENSE.md"><img src="https://img.shields.io/badge/source-BSD--3--Clause-blue.svg?style=flat-square" alt="Source license"></a>
<a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/combined%20work-GPL--3.0-blue.svg?style=flat-square" alt="Combined work: GPL-3.0"></a>
<br>
<a href="https://claude.ai"><img src="https://img.shields.io/badge/Built%20with-Claude%20Code-D97757?logo=anthropic&logoColor=white&style=flat-square" alt="Built with Claude Code"></a>
<a href="https://b612foundation.org/asteroid-institute/"><img src="https://img.shields.io/badge/Asteroid%20Institute-b612foundation.org-1a1a2e?style=flat-square" alt="Asteroid Institute"></a>
<a href="https://dirac.astro.washington.edu/"><img src="https://img.shields.io/badge/DIRAC%20Institute-dirac.astro.washington.edu-1a1a2e?style=flat-square" alt="DIRAC Institute"></a>

---

Implements [`Propagator`](https://github.com/moeyensj/thor-rs-propagator) over
[ASSIST](https://github.com/matthewholman/assist) +
[REBOUND](https://github.com/hannorein/rebound) via
[assist-rs](https://github.com/B612-Asteroid-Institute/assist-rs). This is the
backend a THOR build with **no proprietary components** uses: everything in
the dependency closure is open source, and because the ASSIST/REBOUND stack
is GPL-3.0, any combined binary is governed by the GPL (this crate's own
source is BSD-3-Clause).

Not yet on crates.io (it publishes after `thor-rs-propagator`, its git
dependency). Until then:

```toml
[dependencies]
thor-rs-assist = { git = "https://github.com/moeyensj/thor-rs-assist.git", rev = "<pin>" }
```

## What it does

- **Propagation** — IAS15 (15th-order, adaptive) N-body integration of test
  particles under the Sun, planets, Moon, and 16 massive asteroids from JPL
  DE440 / sb441-n16, through assist-rs. Finite-difference state transition
  matrices drive first-order covariance transport; the integrator defaults
  are tuned to bound runaway step collapse near perihelion (with the
  SIGSEGV-history rationale documented on the tuning constants).
- **Ephemerides** — light-time-iterated astrometric topocentric spherical
  coordinates with rates, analytic observation Jacobians, and the
  STM-at-emit covariance convention THOR's OD consumes.
- **Observatories & data** — MPC observatory handling and ephemeris/kernel
  acquisition through assist-rs's data manager (`~/.cache/assist-rs`):
  `de440.bsp`, `sb441-n16.bsp`, Earth-orientation BPCs, and the extended
  obscodes table, fetched on first use or pre-staged.
- **Capability honesty** — sample-based covariance methods and fast-integrator
  requests it cannot honor are loud typed `Unsupported` errors, per the
  contract crate's rules. `force_model()` reports
  `assist:de440+sb441-n16+ias15-prs23` so every fitted orbit records its
  physics.

## Quick start

```rust,no_run
use thor_rs_propagator::{PropagationConfig, Propagator};

let prop = thor_rs_assist::AssistPropagator::from_data_manager()?;
let states = prop.propagate(&my_test_orbit, &[60810.0], &PropagationConfig::default())?;
# Ok::<(), thor_rs_propagator::PropagatorError>(())
```

`from_data_manager()` downloads the SPICE kernels and obscodes table into
`~/.cache/assist-rs` on first use (network required once). The REBOUND SPK
reader is not thread-safe across concurrent ephemeris loads — construct one
propagator and share it, and run data-touching tests single-threaded.

## Provenance

Extracted 2026-07-30 from THOR v2's `src/propagator/assist.rs` (repository
`moeyensj/thor_rust`), where the adapter and its integrator tuning were
developed and production-tested. THOR consumes this crate behind its
opt-in `assist` feature and keeps its end-to-end assist tests in-tree as a
dev-dependency.

## Acknowledgments

Developed with support from the [Asteroid Institute](https://b612foundation.org/asteroid-institute/)
(a program of the B612 Foundation) and the [DIRAC Institute](https://dirac.astro.washington.edu/)
at the University of Washington.

## License

This crate's source: BSD 3-Clause ([LICENSE.md](LICENSE.md)). Its mandatory
dependencies assist-rs / libassist-sys / librebound-sys are GPL-3.0, so any
distributed binary containing this backend is a GPL-3.0 combined work —
which is exactly why THOR ships it as an opt-in source feature rather than
in prebuilt artifacts.
