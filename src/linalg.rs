//! The small matrix/constant slice the adapter needs from THOR.

/// Radians -> degrees.
pub const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

pub use hyperjet::linalg::mat6::{mat6_mul, mat6_transpose};

fn vec6_norm(x: &[f64; 6]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Find the eigenvector corresponding to the largest-magnitude eigenvalue
/// of a 6×6 symmetric matrix via power iteration.
///
/// Returns (eigenvector, eigenvalue). The eigenvector is unit-length.
/// The eigenvalue is refined via Rayleigh quotient after convergence.
pub fn mat6_eigenvector_max(a: &[[f64; 6]; 6], max_iter: usize, tol: f64) -> ([f64; 6], f64) {
    let inv_sqrt6 = 1.0 / (6.0_f64).sqrt();
    let mut v = [inv_sqrt6; 6];

    #[allow(unused_assignments)]
    let mut lambda = 0.0_f64;

    for _ in 0..max_iter {
        // A v
        let mut w = [0.0; 6];
        for i in 0..6 {
            let mut sum = 0.0;
            for k in 0..6 {
                sum += a[i][k] * v[k];
            }
            w[i] = sum;
        }
        let w_norm = vec6_norm(&w);

        if w_norm < 1e-30 {
            return (v, 0.0);
        }

        let mut dot = 0.0;
        for i in 0..6 {
            dot += w[i] * v[i];
        }
        let sign = if dot >= 0.0 { 1.0 } else { -1.0 };
        let lambda_new = sign * w_norm;

        let mut residual_sq = 0.0;
        for i in 0..6 {
            let r = w[i] - lambda_new * v[i];
            residual_sq += r * r;
        }
        let converged = residual_sq.sqrt() / lambda_new.abs() < tol;

        #[allow(unused_assignments)]
        {
            lambda = lambda_new;
        }

        for i in 0..6 {
            v[i] = w[i] / w_norm;
        }

        if converged {
            break;
        }
    }

    // Rayleigh quotient refinement
    let mut av = [0.0; 6];
    for i in 0..6 {
        let mut sum = 0.0;
        for k in 0..6 {
            sum += a[i][k] * v[k];
        }
        av[i] = sum;
    }
    lambda = 0.0;
    for i in 0..6 {
        lambda += av[i] * v[i];
    }

    (v, lambda)
}

/// Minimum eigenvalue of a symmetric 6×6 matrix.
///
/// Computed as \(c - \lambda_{max}(cI - A)\) with \(c = \lVert A\rVert_F\),
/// which bounds the spectral radius, so \(cI - A\) is positive
/// semidefinite and [`mat6_eigenvector_max`]'s power iteration converges
/// to its dominant eigenvalue. Deterministic, ~exact for the
/// definiteness tagging it serves (covariance quality), not a
/// high-precision spectral routine.
pub fn mat6_min_eigenvalue_sym(a: &[[f64; 6]; 6]) -> f64 {
    let mut frob_sq = 0.0;
    for row in a.iter() {
        for &v in row.iter() {
            frob_sq += v * v;
        }
    }
    let c = frob_sq.sqrt();
    if c == 0.0 {
        return 0.0; // The zero matrix: every eigenvalue is 0.
    }
    let mut shifted = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            shifted[i][j] = -a[i][j];
        }
        shifted[i][i] += c;
    }
    let (_, lambda_max_shifted) = mat6_eigenvector_max(&shifted, 200, 1e-12);
    c - lambda_max_shifted
}

pub fn propagate_covariance_6x6(jacobian: &[[f64; 6]; 6], sigma: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let jt = mat6_transpose(jacobian);
    let j_sigma = mat6_mul(jacobian, sigma);
    mat6_mul(&j_sigma, &jt)
}
