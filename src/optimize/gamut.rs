use super::cie1931::*;

/// Color gamuts
pub enum Gamut {
    SRGB,
    ProPhotoRGB,
    ACES2065_1,
    REC2020,
    ERGB,
    XYZ,
}

pub(crate) struct GamutData {
    pub(crate) illuminant: [f64; SAMPLES],
    pub(crate) xyz_to_rgb: [[f64; 3]; 3],
    pub(crate) rgb_to_xyz: [[f64; 3]; 3],
}

impl Gamut {
    pub fn parse(str: &str) -> Option<Self> {
        if str.eq_ignore_ascii_case("sRGB") {
            Some(Gamut::SRGB)
        } else if str.eq_ignore_ascii_case("eRGB") {
            Some(Gamut::ERGB)
        } else if str.eq_ignore_ascii_case("XYZ") {
            Some(Gamut::XYZ)
        } else if str.eq_ignore_ascii_case("ProPhotoRGB") {
            Some(Gamut::ProPhotoRGB)
        } else if str.eq_ignore_ascii_case("ACES2065_1") {
            Some(Gamut::ACES2065_1)
        } else if str.eq_ignore_ascii_case("REC2020") {
            Some(Gamut::REC2020)
        } else {
            None
        }
    }

    pub(crate) fn get_data(&self) -> GamutData {
        match self {
            Gamut::SRGB => data(ILLUMINANT_D65, XYZ_TO_SRGB, SRGB_TO_XYZ),
            Gamut::ERGB => data(ILLUMINANT_E, XYZ_TO_ERGB, ERGB_TO_XYZ),
            Gamut::XYZ => data(ILLUMINANT_E, XYZ_TO_XYZ, XYZ_TO_XYZ),
            Gamut::ProPhotoRGB => data(ILLUMINANT_D50, XYZ_TO_PROPHOTO_RGB, PROPHOTO_RGB_TO_XYZ),
            Gamut::ACES2065_1 => data(ILLUMINANT_D60, XYZ_TO_ACES2065_1, ACES2065_1_TO_XYZ),
            Gamut::REC2020 => data(ILLUMINANT_D65, XYZ_TO_REC2020, REC2020_TO_XYZ),
        }
    }
}

fn data(
    illuminant: [f64; SAMPLES],
    xyz_to_rgb: [[f64; 3]; 3],
    rgb_to_xyz: [[f64; 3]; 3],
) -> GamutData {
    GamutData {
        illuminant,
        xyz_to_rgb,
        rgb_to_xyz,
    }
}

static XYZ_TO_SRGB: [[f64; 3]; 3] = [
    [3.240479, -1.537150, -0.498535],
    [-0.969256, 1.875991, 0.041556],
    [0.055648, -0.204043, 1.057311],
];

static SRGB_TO_XYZ: [[f64; 3]; 3] = [
    [0.412453, 0.357580, 0.180423],
    [0.212671, 0.715160, 0.072169],
    [0.019334, 0.119193, 0.950227],
];

static XYZ_TO_XYZ: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

static XYZ_TO_ERGB: [[f64; 3]; 3] = [
    [2.689989, -1.276020, -0.413844],
    [-1.022095, 1.978261, 0.043821],
    [0.061203, -0.224411, 1.162859],
];

static ERGB_TO_XYZ: [[f64; 3]; 3] = [
    [0.496859, 0.339094, 0.164047],
    [0.256193, 0.678188, 0.065619],
    [0.023290, 0.113031, 0.863978],
];

static XYZ_TO_PROPHOTO_RGB: [[f64; 3]; 3] = [
    [1.3459433, -0.2556075, -0.0511118],
    [-0.5445989, 1.5081673, 0.0205351],
    [0.0000000, 0.0000000, 1.2118128],
];

static PROPHOTO_RGB_TO_XYZ: [[f64; 3]; 3] = [
    [0.7976749, 0.1351917, 0.0313534],
    [0.2880402, 0.7118741, 0.0000857],
    [0.0000000, 0.0000000, 0.8252100],
];

static XYZ_TO_ACES2065_1: [[f64; 3]; 3] = [
    [1.0498110175, 0.0000000000, -0.0000974845],
    [-0.4959030231, 1.3733130458, 0.0982400361],
    [0.0000000000, 0.0000000000, 0.9912520182],
];

static ACES2065_1_TO_XYZ: [[f64; 3]; 3] = [
    [0.9525523959, 0.0000000000, 0.0000936786],
    [0.3439664498, 0.7281660966, -0.0721325464],
    [0.0000000000, 0.0000000000, 1.0088251844],
];

static XYZ_TO_REC2020: [[f64; 3]; 3] = [
    [1.7166511880, -0.3556707838, -0.2533662814],
    [-0.6666843518, 1.6164812366, 0.0157685458],
    [0.0176398574, -0.0427706133, 0.9421031212],
];

static REC2020_TO_XYZ: [[f64; 3]; 3] = [
    [0.6369580483, 0.1446169036, 0.1688809752],
    [0.2627002120, 0.6779980715, 0.0593017165],
    [0.0000000000, 0.0280726930, 1.0609850577],
];
