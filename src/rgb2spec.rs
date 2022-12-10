//! # rgb2spec-rs
//!
//! This crate is a port of [rgb2spec](https://github.com/mitsuba-renderer/rgb2spec),
//! an implementation of [this paper](http://rgl.epfl.ch/publications/Jakob2019Spectral)
//! by Wenzel Jakob and Johannes Hanika.
//!
//! It can be used to convert RGB colors in various color spaces to coefficent representations of reflectance spectra.
//! These spectra can then be evaluated for wavelengths in the visible spectrum.
//!
//! # Example
//!
//! ```rust
//! use rgb2spec::{self, RGB2Spec};
//!
//! let rgb = [1.0, 1.0, 0.3];
//! let rgb2spec = RGB2Spec::load("examples/out.spec")?;
//!
//! let coefficients = rgb2spec.fetch(rgb);
//! let wavelength = 480.0;
//! let reflectance = rgb2spec::eval_precise(coefficients, wavelength);
//! ```

use std::fmt;
use std::fs::File;
use std::io::{Error, ErrorKind, Read, Write};
use std::path::Path;

pub mod optimize;

/// Start of the wavelength range that this crate has data for.
pub const LAMBDA_MIN: f64 = optimize::cie1931::LAMBDA_MIN;
/// End of the wavelength range that this crate has data for.
pub const LAMBDA_MAX: f64 = optimize::cie1931::LAMBDA_MAX;
/// Size of the wavelength range that this crate has data for.
pub const LAMBDA_RANGE: f64 = optimize::cie1931::LAMBDA_RANGE;

const N_COEFFS: usize = 3;

/// A precomputed model used to convert RGB data to a coefficient representation of reflectance spectra.
///
/// This crate provides the following methods of instantiating this struct:
/// * Using [RGB2Spec::load] to load the model from a file
/// * Using [optimize](optimize::optimize) to compute the model (slow)
///
/// The crate also includes a CLI program that can be used to compute a model and save it to a file.
/// Use `cargo run` in the crate's root to execute it.
pub struct RGB2Spec {
    resolution: u32,
    scale: Vec<f32>,
    data: Vec<f32>,
}

impl RGB2Spec {
    /// Loads a [RGB2Spec] model from a file.
    ///
    /// The binary format is compatible with the [original implementation](https://github.com/mitsuba-renderer/rgb2spec).
    ///
    /// Because the binary format doesn't contain information on which [Gamut](optimize::Gamut)
    /// was used to generate the model there is no interface to retrieve this.
    /// When this information is required the user will need to keep track of it manually.
    ///
    /// Returns a [std::io::Error] if the file cannot be opened or does not comply with the format.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<RGB2Spec, Error> {
        let mut file = File::open(path)?;
        RGB2Spec::from_reader(&mut file)
    }

    /// Loads a [RGB2Spec] model from a reader.
    ///
    /// The binary format is compatible with the [original implementation](https://github.com/mitsuba-renderer/rgb2spec).
    ///
    /// Because the binary format doesn't contain information on which [Gamut](optimize::Gamut)
    /// was used to generate the model there is no interface to retrieve this.
    /// When this information is required the user will need to keep track of it manually.
    ///
    /// Returns a [std::io::Error] if the reader cannot be read or does not comply with the format.
    pub fn from_reader<R: Read>(reader: &mut R) -> Result<RGB2Spec, Error> {
        let mut header = [0; 4];
        reader.read_exact(&mut header)?;

        if std::str::from_utf8(&header) != Ok("SPEC") {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "Header is not correct. (Expected SPEC as bytes)",
            ));
        }

        let mut resolution = [0; 4];
        reader.read_exact(&mut resolution)?;
        let resolution = u32::from_le_bytes(resolution);

        let res = resolution as usize;

        // Prevent overflows and OOMs
        if res > 0xFFFF {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "Parsed resolution too large.",
            ));
        }

        let data_len = res * res * res * 3 * N_COEFFS;

        let mut scale = vec![0_u8; res * std::mem::size_of::<f32>()];
        let mut data = vec![0_u8; data_len * std::mem::size_of::<f32>()];

        reader.read_exact(&mut scale)?;
        reader.read_exact(&mut data)?;

        let scale = scale
            .chunks_exact(std::mem::size_of::<f32>())
            .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect::<Vec<_>>();

        assert_eq!(scale.len(), res);

        let data = data
            .chunks_exact(std::mem::size_of::<f32>())
            .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect::<Vec<_>>();

        assert_eq!(data.len(), data_len);

        Ok(RGB2Spec {
            resolution,
            scale,
            data,
        })
    }

    /// Saves the model to a file.
    ///
    /// The binary format is compatible with the [original implementation](https://github.com/mitsuba-renderer/rgb2spec).
    ///
    /// Returns a [std::io::Error] if the file cannot be written to.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        let mut file = File::create(path)?;
        self.to_writer(&mut file)
    }

    /// Writes the model to a writer.
    ///
    /// The binary format is compatible with the [original implementation](https://github.com/mitsuba-renderer/rgb2spec).
    ///
    /// Returns a [std::io::Error] if the writer cannot be written to.
    pub fn to_writer<W: Write>(&self, writer: &mut W) -> Result<(), Error> {
        writer.write_all("SPEC".as_bytes())?;
        writer.write_all(&self.resolution.to_le_bytes())?;

        writer.write_all(
            &self
                .scale
                .iter()
                .flat_map(|&x| x.to_le_bytes().into_iter())
                .collect::<Vec<_>>(),
        )?;

        writer.write_all(
            &self
                .data
                .iter()
                .flat_map(|&x| x.to_le_bytes().into_iter())
                .collect::<Vec<_>>(),
        )?;

        Ok(())
    }

    /// Convert an RGB tuple into a RGB2Spec coefficient representation.
    ///
    /// The spectrum that the coefficients represent can then be sampled by passing the coefficents to [eval_precise].
    pub fn fetch(&self, rgb: [f32; 3]) -> [f32; 3] {
        // Determine largest RGB component
        let res = self.resolution;
        let rgb = rgb.map(|x| x.clamp(0.0, 1.0));

        let mut i = 0;

        for j in 1..3 {
            if rgb[j] >= rgb[i] {
                i = j;
            }
        }

        let z = rgb[i];
        // Prevent NaN values for (0, 0, 0)
        let scale = if z > 0.0 { (res - 1) as f32 / z } else { 0.0 };
        let x = rgb[(i + 1) % 3] * scale;
        let y = rgb[(i + 2) % 3] * scale;

        // Trilinearly interpolated lookup
        let xi = (x as u32).min(res - 2);
        let yi = (y as u32).min(res - 2);
        let zi = find_interval(&self.scale, self.resolution, z);
        let mut offset = (((i as u32 * res + zi) * res + yi) * res + xi) as usize * N_COEFFS;
        let dx = N_COEFFS;
        let dy = (N_COEFFS as u32 * res) as usize;
        let dz = (N_COEFFS as u32 * res * res) as usize;

        let x1 = x - xi as f32;
        let x0 = 1.0 - x1;
        let y1 = y - yi as f32;
        let y0 = 1.0 - y1;
        let z1 =
            (z - self.scale[zi as usize]) / (self.scale[zi as usize + 1] - self.scale[zi as usize]);
        let z0 = 1.0 - z1;

        let mut out = [0.0; N_COEFFS];

        for o in &mut out {
            *o = ((self.data[offset] * x0 + self.data[offset + dx] * x1) * y0
                + (self.data[offset + dy] * x0 + self.data[offset + dy + dx] * x1) * y1)
                * z0
                + ((self.data[offset + dz] * x0 + self.data[offset + dz + dx] * x1) * y0
                    + (self.data[offset + dz + dy] * x0 + self.data[offset + dz + dy + dx] * x1)
                        * y1)
                    * z1;
            offset += 1;
        }

        out
    }

    pub(crate) fn new(resolution: u32, scale: Vec<f32>, data: Vec<f32>) -> Self {
        RGB2Spec {
            resolution,
            scale,
            data,
        }
    }
}

impl fmt::Debug for RGB2Spec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RGB2Spec {{ resolution: {}, scale: Vec<f32> {{ len: {} }}, data: Vec<f32> {{ len: {} }} }}",
            self.resolution,
            self.scale.len(),
            self.data.len(),
        )
    }
}

/// Used to evaluate the spectrum represented by coefficients obtained from [RGB2Spec::fetch] for a specific wavelength.
///
/// Results are only accurate for wavelengths in the range [[LAMBDA_MIN], [LAMBDA_MAX]].
#[inline]
pub fn eval_precise(coefficients: [f32; 3], wavelength: f32) -> f32 {
    let x = fma(
        fma(coefficients[0], wavelength, coefficients[1]),
        wavelength,
        coefficients[2],
    );
    let y = 1.0 / fma(x, x, 1.0).sqrt();
    fma(0.5 * x, y, 0.5)
}

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a * b + c
}

fn find_interval(values: &[f32], size: u32, x: f32) -> u32 {
    let mut left = 0;
    let last_interval = size - 2;
    let mut size = last_interval;

    while size > 0 {
        let half = size >> 1;
        let middle = left + half + 1;

        if values[middle as usize] <= x {
            left = middle;
            size -= half + 1;
        } else {
            size = half;
        }
    }

    left.min(last_interval)
}
