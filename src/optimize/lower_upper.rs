#![allow(clippy::needless_range_loop)]

use std::error::Error;

const N: usize = 3;

// LU decomposition & triangular solving code lifted from Wikipedia

/* INPUT: a - square matrix to decompose
 *        tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Decomposed matrix, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */
pub(crate) fn decompose(
    a: &[[f64; 3]; 3],
    tol: f64,
) -> Result<([[f64; 3]; 3], [usize; 4]), DegenerateMatrixError> {
    let mut a = *a;

    //Unit permutation matrix, p[N] initialized with N
    let mut p = [0, 1, 2, 3];

    for i in 0..N {
        let mut max_a = 0.0;
        let mut i_max = i;

        for k in i..N {
            let abs_a = a[k][i].abs();

            if abs_a > max_a {
                max_a = abs_a;
                i_max = k;
            }
        }

        if max_a < tol {
            return Err(DegenerateMatrixError {});
        }

        if i_max != i {
            p.swap(i, i_max);
            a.swap(i, i_max);

            //counting pivots starting from N (for determinant)
            p[N] += 1;
        }

        for j in (i + 1)..N {
            a[j][i] /= a[i][i];

            for k in (i + 1)..N {
                a[j][k] -= a[j][i] * a[i][k];
            }
        }
    }

    Ok((a, p))
}

#[derive(Debug)]
pub struct DegenerateMatrixError {}

impl std::fmt::Display for DegenerateMatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Decomposition failed, matrix is degenerate.")
    }
}

impl Error for DegenerateMatrixError {}

/* INPUT: a,p filled in decompose; b - rhs vector
 * OUTPUT: x - solution vector of a*x=b
 */
pub(crate) fn solve(a: [[f64; 3]; 3], p: [usize; 4], b: [f64; 3]) -> [f64; 3] {
    let mut x = [0.0; 3];

    for i in 0..N {
        x[i] = b[p[i]];

        for k in 0..i {
            x[i] -= a[i][k] * x[k];
        }
    }

    for i in (0..=(N - 1)).rev() {
        for k in (i + 1)..N {
            x[i] -= a[i][k] * x[k];
        }

        x[i] /= a[i][i];
    }

    x
}
