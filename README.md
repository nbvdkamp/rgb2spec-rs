# rgb2spec-rs

This crate is a port of [rgb2spec](https://github.com/mitsuba-renderer/rgb2spec), an implementation of the paper ["A Low-Dimensional Function Space for Efficient Spectral Upsampling"](http://rgl.epfl.ch/publications/Jakob2019Spectral) by Wenzel Jakob and Johannes Hanika.

It can be used to convert RGB colors in various color spaces to coefficent representations of reflectance spectra. These spectra can then be evaluated for wavelengths in the visible spectrum.

### Usage

See the [crate documentation](https://docs.rs/) or `examples/` for example usages.

### CLI
The crate also includes a command line program that can be used to compute a model and save it to a file. Use `cargo run` in the crate's root to execute it.