//! HEALPix projection.

use std::f64::consts::PI;

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};
use crate::math::HALF_PI;

/// Mask to keep only the f64 sign
pub const F64_SIGN_BIT_MASK: u64 = 0x8000000000000000;
/// Equals !F64_SIGN_BIT_MASK (the inverse of the f64 sign mask)
pub const F64_BUT_SIGN_BIT_MASK: u64 = 0x7FFFFFFFFFFFFFFF;

/// Limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
/// Equals 2/3, see Eq. (1) in Gorsky2005.
pub const TRANSITION_Z: f64 = 2_f64 / 3_f64;
/// Inverse of the limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
/// Equals 1/(2/3) = 1.5, see Eq. (1) in Gorsky2005.
pub const ONE_OVER_TRANSITION_Z: f64 = 1.5_f64;

/// Upper limit on sqrt(3(1-|z|)) to consider that we are not near from the poles
const EPS_POLE: f64 = 1e-13_f64;

const ONE_OVER_SQRT6: f64 = 0.408_248_290_463_863_f64;

pub const PI_OVER_FOUR: f64 = PI / 4.0;
pub const FOUR_OVER_PI: f64 = 4.0 / PI;


/// HEALPix projection.
pub struct Hpx;

impl Default for Hpx {
  fn default() -> Self {
    Self::new()
  }
}

impl Hpx {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Hpx {
  
  const NAME: &'static str = "HEALPix";
  const WCS_NAME: &'static str = "HPX";
  
  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-PI..=PI),
      Some(-HALF_PI..=HALF_PI)
    );
    &PROJ_BOUNDS
  }

  /// Returns the projection, in the 2D Euclidean plane, of the given position on the unit sphere.
  /// # Output
  /// * `X`: coordinate along the X-axis in the projection plane, in `[-PI, PI]`
  /// * `Y`: coordinate along the Y-axis in the projection plane, in `[-PI/2, PI/2]`
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if xyz.z > TRANSITION_Z {
      // North polar cap, Collignon projection.
      let (x_pm1, offset) = xpm1_and_offset(xyz.x, xyz.y);
      let sqrt_3_one_min_z = (3.0 * one_minus_z_pos(xyz)).sqrt();
      Some(ProjXY::new(
        ((x_pm1 * sqrt_3_one_min_z) + offset as f64) * PI_OVER_FOUR,
        (2.0 - sqrt_3_one_min_z) * PI_OVER_FOUR
      ))
    } else if xyz.z < -TRANSITION_Z {
      // South polar cap, Collignon projection
      let (x_pm1, offset) = xpm1_and_offset(xyz.x, xyz.y);
      let sqrt_3_one_min_z = (3.0 * one_minus_z_neg(xyz)).sqrt();
      Some(ProjXY::new(
        ((x_pm1 * sqrt_3_one_min_z) + offset as f64) * PI_OVER_FOUR,
        (-2.0 + sqrt_3_one_min_z) * PI_OVER_FOUR
      ))
    } else {
      // Equatorial region, Cylindrical equal area projection
      Some(ProjXY::new(
        xyz.y.atan2(xyz.x),
        xyz.z * ONE_OVER_TRANSITION_Z * PI_OVER_FOUR
      ))
    }
  }

  /// # Input
  /// * `X`: coordinate along the X-axis in the projection plane, in `[-PI, PI]`
  /// * `Y`: coordinate along the Y-axis in the projection plane, in `[-PI/2, PI/2]`
  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let mut z = pos.y * FOUR_OVER_PI;
    let x = pos.x * FOUR_OVER_PI;
    if !(-2f64..=2f64).contains(&z) || !(-4f64..4f64).contains(&x) {
      None
    } else if z > 1.0 {
      // North polar cap
      let x = abs_sign_decompose(x);
      let OffsetAndPM1 { offset, mut pm1 } = pm1_offset_decompose(x.abs);
      deproj_collignon(&mut pm1, &mut z);
      if (-1.0..=1.0).contains(&pm1) {
        apply_offset_and_signs(&mut pm1, offset, x.sign);
        pm1 *= PI_OVER_FOUR;
        let (sinb, cosb) = z.sin_cos();
        let (sinl, cosl) = pm1.sin_cos();
        Some(XYZ::new(cosl * cosb, sinl * cosb, sinb))
      } else {
        None
      }
    } else if z < -1.0 {
      // South polar cap
      let x = abs_sign_decompose(x);
      let OffsetAndPM1 { offset, mut pm1 } = pm1_offset_decompose(x.abs);
      z = -z;
      deproj_collignon(&mut pm1, &mut z);
      if (-1.0..=1.0).contains(&pm1) {
        apply_offset_and_signs(&mut pm1, offset, x.sign);
        pm1 *= PI_OVER_FOUR;
        let (sinb, cosb) = (-z).sin_cos();
        let (sinl, cosl) = pm1.sin_cos();
        Some(XYZ::new(cosl * cosb, sinl * cosb, sinb))
      } else {
        None
      }
    }  else {
      // Equatorial region
      let z = z * TRANSITION_Z; // z = sin(lat) = sinb
      let cosb = (1.0 - z.pow2()).sqrt();
      let (sinl, cosl) = pos.x.sin_cos();
      Some(XYZ::new(cosl * cosb, sinl * cosb, z))
    }
  }
}

fn xpm1_and_offset(x: f64, y: f64) -> (f64, i8) {
  let x_neg = (x < 0.0) as i8;             debug_assert!(x_neg == 0 || x_neg == 1);
  let y_neg = (y < 0.0) as i8;             debug_assert!(y_neg == 0 || y_neg == 1);
  // x>0, y>0 => [    0,  pi/2[ => offset =  1
  // x<0, y>0 => [pi/2 ,    pi[ => offset =  3
  // x<0, y<0 => [3pi/2,    pi[ => offset = -3
  // x>0, y<0 => [pi   , 3pi/2[ => offset = -1
  let offset = ((-y_neg) << 2) + 1 + ((x_neg ^ y_neg) << 1);
  let lon = y.abs().atan2(x.abs()); debug_assert!((0.0..=PI / 2.0).contains(&lon));
  let x02 = lon * FOUR_OVER_PI;           debug_assert!((0.0..=2.0).contains(&x02));
  if x_neg != y_neg { // Could be replaced by a sign copy from (x_neg ^ y_neg) << 32
    (1.0 - x02, offset)
  } else {
    (x02 - 1.0, offset)
  }
}

fn one_minus_z_pos(xyz: &XYZ) -> f64 {
  debug_assert!(xyz.z > 0.0);
  let d2 = xyz.x.pow2() + xyz.y.pow2(); // z = sqrt(1 - d2) AND sqrt(1 - x) = 1 - x / 2 - x^2 / 8 - x^3 / 16 - 5 x^4/128 - 7 * x^5/256
  if d2 < 1.0e-3 { // <=> dec > 88.187846253 deg
    d2 * (0.5 + d2 * (0.125 + d2 * (0.0625 + d2 * (0.0390625 + d2 * 0.02734375))))
  } else {
    1.0 - xyz.z
  }
}

fn one_minus_z_neg(xyz: &XYZ) -> f64 {
  debug_assert!(xyz.z < 0.0);
  let d2 = xyz.x.pow2() + xyz.y.pow2(); // z = sqrt(1 - d2) AND sqrt(1 - x) = 1 - x / 2 - x^2 / 8 - x^3 / 16 - 5 x^4/128 - 7 * x^5/256
  if d2 < 1.0e-3 { // <=> dec < -88.187846253 deg
    // 0.5 * d2 + 0.125 * d2 * d2
    d2 * (0.5 + d2 * (0.125 + d2 * (0.0625 + d2 * (0.0390625 + d2 * 0.02734375))))
  } else {
    xyz.z + 1.0
  }
}

// Returns the absolute value of the given double together with its bit of sign
struct AbsAndSign {
  abs: f64,
  sign: u64,
}
fn abs_sign_decompose(x: f64) -> AbsAndSign {
  let bits = f64::to_bits(x);
  AbsAndSign {
    abs: f64::from_bits(bits & F64_BUT_SIGN_BIT_MASK),
    sign: bits & F64_SIGN_BIT_MASK,
  }
}

// Decompose the given positive real value in
// --* an integer offset in [1, 3, 5, 7] (*PI/4) and
// --* a real value in [-1.0, 1.0] (*PI/4)
struct OffsetAndPM1 {
  offset: u8, // = 1, 3, 5 or 7
  pm1: f64,   // in [-1.0, 1.0]
}
fn pm1_offset_decompose(x: f64) -> OffsetAndPM1 {
  let floor: u8 = x as u8;
  let odd_floor: u8 = floor | 1u8;
  OffsetAndPM1 {
    offset: odd_floor & 7u8, // value modulo 8
    pm1: x - (odd_floor as f64),
  }
}

fn deproj_collignon(lon: &mut f64, lat: &mut f64) {
  *lat = 2.0 - *lat;
  if is_not_near_from_pole(*lat) { // Rare, so few risks of branch miss-prediction
    *lon /= *lat;
    deal_with_numerical_approx_in_edges(lon);
  } // in case of pole, lon = lat = 0 (we avoid NaN due to division by lat=0)
  *lat *= ONE_OVER_SQRT6;
  // Using acos is OK here since lat < 1/sqrt(6), so not near from 1.
  *lat = 2.0 * (f64::acos(*lat) - PI_OVER_FOUR);
}

fn is_not_near_from_pole(sqrt_of_three_time_one_minus_sin_of: f64) -> bool {
  // In case of pole: x = y = 0
  sqrt_of_three_time_one_minus_sin_of > EPS_POLE
}

fn deal_with_numerical_approx_in_edges(lon: &mut f64) {
  if *lon > 1.0 {
    *lon = 1.0;
  } else if *lon < -1.0 {
    *lon = -1.0;
  }
}

// Shift x by the given offset and apply lon sign to x
pub(crate) fn apply_offset_and_signs(lon: &mut f64, off: u8, lon_sign: u64) {
  *lon += off as f64;
  *lon = f64::from_bits(f64::to_bits(*lon) | lon_sign);
}