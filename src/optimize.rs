#![allow(clippy::needless_range_loop)]

use std::error::Error;
use std::fmt;

pub(crate) mod cie1931;
pub mod gamut;
mod lower_upper;

use super::RGB2Spec;
use cie1931 as CIE;
use gamut::*;

/// Discretization of quadrature scheme
const CIE_FINE_SAMPLES: usize = (CIE::SAMPLES - 1) * 3 + 1;
const EPSILON: f64 = 1e-4;

/// Computes an [RGB2Spec] model for the given resolution and [Gamut].
///
/// Returns an error when the Gauss-Newton algorithm for optimizing the coefficients runs into a degenerate matrix.
/// This occurs for low resolutions (less than 10).
pub fn optimize(gamut: Gamut, resolution: usize) -> Result<RGB2Spec, MathError> {
    let tables = init_tables(gamut);

    let scale = (0..resolution)
        .map(|k| smoothstep(smoothstep(k as f64 / (resolution - 1) as f64)) as f32)
        .collect::<Vec<_>>();

    let buffer_size = super::N_COEFFS * 3 * resolution * resolution * resolution;
    let mut data = vec![0.0; buffer_size];

    for l in 0..3 {
        for j in 0..resolution {
            let y = j as f64 / (resolution - 1) as f64;

            for i in 0..resolution {
                let x = i as f64 / (resolution - 1) as f64;
                let mut rgb = [0.0; 3];
                let start = resolution / 5;

                let mut iterate = |start, end| {
                    let mut coeffs = [0.0; 3];

                    // Manual for loop because range and reverse range types aren't compatible
                    let increment = if start > end { -1 } else { 1 };
                    let mut k = start;

                    while k != end {
                        let b = scale[k as usize] as f64;

                        rgb[l] = b;
                        rgb[(l + 1) % 3] = x * b;
                        rgb[(l + 2) % 3] = y * b;

                        coeffs = gauss_newton(rgb, coeffs, &tables)?;

                        let c0 = CIE::LAMBDA_MIN;
                        let c1 = 1.0 / CIE::LAMBDA_RANGE;
                        let a = coeffs[0];
                        let b = coeffs[1];
                        let c = coeffs[2];

                        let idx = ((l * resolution + k as usize) * resolution + j) * resolution + i;

                        data[3 * idx] = (a * (sqr(c1))) as f32;
                        data[3 * idx + 1] = (b * c1 - 2.0 * a * c0 * (sqr(c1))) as f32;
                        data[3 * idx + 2] = (c - b * c0 * c1 + a * (sqr(c0 * c1))) as f32;

                        k += increment;
                    }
                    Ok(())
                };

                iterate(start as i64, resolution as i64)?;
                iterate(start as i64, -1)?;
            }
        }
    }

    Ok(RGB2Spec::new(resolution as u32, scale, data))
}

#[derive(Debug)]
pub struct MathError {
    source: lower_upper::DegenerateMatrixError,
    rgb: [f64; 3],
    coeffs: [f64; 3],
}

impl fmt::Display for MathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LU decomposition failed! rgb: {:?}, coeffs: {:?}",
            self.rgb, self.coeffs
        )
    }
}

impl Error for MathError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}

fn sigmoid(x: f64) -> f64 {
    0.5 * x / (1.0 + x * x).sqrt() + 0.5
}

fn smoothstep(x: f64) -> f64 {
    x * x * (3.0 - 2.0 * x)
}

fn sqr(x: f64) -> f64 {
    x * x
}

fn cie_lab(p: [f64; 3], data: &Data) -> [f64; 3] {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    let xw = data.xyz_whitepoint[0];
    let yw = data.xyz_whitepoint[1];
    let zw = data.xyz_whitepoint[2];

    for j in 0..3 {
        x += p[j] * data.rgb_to_xyz[0][j];
        y += p[j] * data.rgb_to_xyz[1][j];
        z += p[j] * data.rgb_to_xyz[2][j];
    }

    let f = |t: f64| {
        let delta: f64 = 6.0 / 29.0;

        if t > delta * delta * delta {
            t.cbrt()
        } else {
            t / (delta * delta * 3.0) + (4.0 / 29.0)
        }
    };

    [
        116.0 * f(y / yw) - 16.0,
        500.0 * (f(x / xw) - f(y / yw)),
        200.0 * (f(y / yw) - f(z / zw)),
    ]
}

/// Precomputed tables for fast spectral -> RGB conversion
struct Data {
    rgb_to_xyz: [[f64; 3]; 3],
    lambda_tbl: [f64; CIE_FINE_SAMPLES],
    rgb_tbl: [[f64; CIE_FINE_SAMPLES]; 3],
    xyz_whitepoint: [f64; 3],
}

///  This function precomputes tables used to convert arbitrary spectra
///  to RGB
///
///  A composite quadrature rule integrates the CIE curves, reflectance, and
///  illuminant spectrum over each 5nm segment in the 360..830nm range using
///  Simpson's 3/8 rule (4th-order accurate), which evaluates the integrand at
///  four positions per segment. While the CIE curves and illuminant spectrum are
///  linear over the segment, the reflectance could have arbitrary behavior,
///  hence the extra precautions.
fn init_tables(gamut: Gamut) -> Data {
    let gamut_data = gamut.get_data();

    let mut lambda_tbl = [0.0; CIE_FINE_SAMPLES];
    let mut rgb_tbl = [[0.0; CIE_FINE_SAMPLES]; 3];
    let mut xyz_whitepoint = [0.0; 3];

    let h = CIE::LAMBDA_RANGE / (CIE_FINE_SAMPLES - 1) as f64;

    for i in 0..CIE_FINE_SAMPLES {
        let lambda = CIE::LAMBDA_MIN + i as f64 * h;

        let xyz = [
            CIE::cie_interp(CIE::OBSERVER_X, lambda),
            CIE::cie_interp(CIE::OBSERVER_Y, lambda),
            CIE::cie_interp(CIE::OBSERVER_Z, lambda),
        ];

        let illuminance = CIE::cie_interp(gamut_data.illuminant, lambda);

        let mut weight = 3.0 / 8.0 * h;

        weight *= if i == 0 || i == CIE_FINE_SAMPLES - 1 {
            1.0
        } else if (i - 1) % 3 == 2 {
            2.0
        } else {
            3.0
        };

        lambda_tbl[i] = lambda;

        for k in 0..3 {
            for j in 0..3 {
                rgb_tbl[k][i] += gamut_data.xyz_to_rgb[k][j] * xyz[j] * illuminance * weight;
            }
        }

        for i in 0..3 {
            xyz_whitepoint[i] += xyz[i] * illuminance * weight;
        }
    }

    Data {
        rgb_to_xyz: gamut_data.rgb_to_xyz,
        lambda_tbl,
        rgb_tbl,
        xyz_whitepoint,
    }
}

fn eval_residual(coeffs: [f64; 3], rgb: [f64; 3], data: &Data) -> [f64; 3] {
    let mut out = [0.0; 3];

    for i in 0..CIE_FINE_SAMPLES {
        // Scale lambda to 0..1 range
        let lambda = (data.lambda_tbl[i] - CIE::LAMBDA_MIN) / CIE::LAMBDA_RANGE;

        // Polynomial
        let mut x = 0.0;

        for c in coeffs {
            x = x * lambda + c;
        }

        // Integrate against precomputed curves
        for j in 0..3 {
            out[j] += data.rgb_tbl[j][i] * sigmoid(x);
        }
    }

    let out = cie_lab(out, data);
    let residual = cie_lab(rgb, data);

    [
        residual[0] - out[0],
        residual[1] - out[1],
        residual[2] - out[2],
    ]
}

fn eval_jacobian(coeffs: [f64; 3], rgb: [f64; 3], data: &Data) -> [[f64; 3]; 3] {
    let mut jac = [[0.0; 3]; 3];

    for i in 0..3 {
        let mut tmp = coeffs;
        tmp[i] -= EPSILON;
        let r0 = eval_residual(tmp, rgb, data);

        let mut tmp = coeffs;
        tmp[i] += EPSILON;
        let r1 = eval_residual(tmp, rgb, data);

        for j in 0..3 {
            jac[j][i] = (r1[j] - r0[j]) * (1.0 / (2.0 * EPSILON));
        }
    }

    jac
}

fn gauss_newton(rgb: [f64; 3], coeffs: [f64; 3], data: &Data) -> Result<[f64; 3], MathError> {
    let mut r;
    let mut coeffs = coeffs;
    let it = 15;

    for _ in 0..it {
        let residual = eval_residual(coeffs, rgb, data);
        let jacobian = eval_jacobian(coeffs, rgb, data);

        let (jacobian, permutation) = match lower_upper::decompose(&jacobian, 1e-15) {
            Ok(t) => t,
            Err(source) => {
                return Err(MathError {
                    rgb,
                    coeffs,
                    source,
                });
            }
        };

        let x = lower_upper::solve(jacobian, permutation, residual);

        r = 0.0;

        for j in 0..3 {
            coeffs[j] -= x[j];
            r += residual[j] * residual[j];
        }

        let max = coeffs[0].max(coeffs[1]).max(coeffs[2]);

        if max > 200.0 {
            for c in &mut coeffs {
                *c *= 200.0 / max;
            }
        }

        if r < 1e-6 {
            break;
        }
    }

    Ok(coeffs)
}
