use rgb2spec::{RGB2Spec, LAMBDA_MIN, LAMBDA_RANGE};

pub const SAMPLES: usize = 95;

/// Load an [RGB2Spec] model from disk and use it to print a table of wavelengths and reflectance values
fn main() {
    match RGB2Spec::load("examples/out.spec") {
        Ok(rgb2spec) => {
            let rgb = [1.0, 0.0, 0.0];
            let coefficients = rgb2spec.fetch(rgb);
            println!("Coefficients: {coefficients:?}");

            for i in 0..SAMPLES {
                let lambda = LAMBDA_MIN + i as f64 / (SAMPLES - 1) as f64 * LAMBDA_RANGE;
                println!(
                    "{lambda},{}",
                    rgb2spec::eval_precise(coefficients, lambda as f32)
                );
            }
        }
        Err(e) => {
            println!("Something went wrong: {}", e);
        }
    }
}
