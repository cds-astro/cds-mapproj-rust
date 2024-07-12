
use std::{
  f64::consts::PI,
  ops::RangeInclusive
};

pub mod math;
pub mod sip;
pub mod img2proj;
pub mod img2celestial;

pub mod zenithal;
pub mod cylindrical;
pub mod pseudocyl;
pub mod conic;
// pub mod pconic; // pseudoconic + polyconic ??
// pub mod quadcube; ??
pub mod hybrid;

pub use math::CustomFloat;

/// Equatorial coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct LonLat {
  lon: f64,
  lat: f64,  
}
impl LonLat {
  
  /// New struct from lon/lat coordinates (in radians).
  /// # Warning
  /// * No tes perrformed so far to ensure that
  ///     + `lon` in `[0, 2pi[`
  ///     + `lat` in `[-pi/2, pi/2]`
  pub fn new(lon: f64, lat: f64) -> Self {
    // TODO: perform checks here!
    Self {lon, lat}
  }
  
  /// Get the longitude
  pub fn lon(&self) -> f64 {
    self.lon
  }

  /// Get the latitude
  pub fn lat(&self) -> f64 {
    self.lat
  }
  
  /// Transform into Euclidean coordinates 
  pub fn to_xyz(&self) -> XYZ {
    let (sinl, cosl) = self.lon.sin_cos();
    let (sinb, cosb) = self.lat.sin_cos();
    XYZ::new(cosl * cosb, sinl * cosb, sinb)
  }

  /// Compute the `Haversine` angular distance (adapterd for small distances but not for very large distance).
  pub fn haversine_dist(&self, rhs: &LonLat) -> f64 {
    let shs = squared_half_segment(
      rhs.lon - self.lon, rhs.lat - self.lat,
      self.lat.cos(), rhs.lat.cos());
    sphe_dist(shs)
  }
}

/// Returns the angular distance corresponding to the given squared half great-circle arc segment
fn sphe_dist(squared_half_segment: f64) -> f64 {
  squared_half_segment.sqrt().asin().twice()
}

/// Returns `(s/2)^2` with `s` the segment (i.e. the Euclidean distance) between 
/// the two given points  `P1` and `P2` on the unit-sphere.
/// We recall that `s = 2 sin(ad/2)` with `ad` the angular distance between the two points.
/// # Input
/// - `dlon` the longitude difference, i.e. (P2.lon - P1.lon), in radians
/// - `dlat` the latitude difference, i.e. (P2.lat - P1.lat), in radians
/// - `cos_lat1` cosine of the latitude of the first point
/// - `cos_lat2` cosine of the latitude of the second point
fn squared_half_segment(dlon: f64, dlat: f64, cos_lat1: f64, cos_lat2: f64) -> f64 {
  dlat.half().sin().pow2() + cos_lat1 * cos_lat2 * dlon.half().sin().pow2()
}

/// Euclidean coordinates (supposedly normalized to 1).
#[derive(Debug, Clone, PartialEq)]
pub struct XYZ {
  x: f64,
  y: f64,
  z: f64,
}
impl XYZ {
  
  /// We assume the norm of the input vector is 1.
  pub fn new(x: f64, y: f64, z: f64) -> Self {
    debug_assert!((-1.0..=1.0).contains(&x), "x: {}; y: {}; z: {}", x, y, z);
    debug_assert!((-1.0..=1.0).contains(&y), "y: {}; y: {}; z: {}", x, y, z);
    debug_assert!((-1.0..=1.0).contains(&z), "x: {}; y: {}; z: {}", x, y, z);
    debug_assert!((1.0 - (x.pow2() + y.pow2() + z.pow2())).abs() < 1e-15);
    Self { x, y, z }
  }

  /// Renormalize the input parameters to ensure the norm of the vector equals 1.
  pub fn new_renorming_if_necessary(x: f64, y: f64, z: f64) -> Self {
    let n = (x.pow2() + y.pow2() + z.pow2()).sqrt();
   if !(0.99999999999999..=1.0).contains(&n) {
     Self { x: x / n, y: y / n, z: z / n }
   } else {
     Self { x, y, z } 
   }
  }

  /// Get the x coordinate
  pub fn x(&self) -> f64 {
    self.x
  }

  /// Get the y coordinate
  pub fn y(&self) -> f64 {
    self.y
  }

  /// Get the z coordinate
  pub fn z(&self) -> f64 {
    self.z
  }
  
  /// Transform into equatorial coordinates
  fn to_lonlat(&self) -> LonLat {
    let r2 = self.x.pow2() + self.y.pow2();
    // Latitude in [-pi/2, pi/2] (ok, since cos always positive here)
    let lat = self.z.atan2(r2.sqrt());
    // Compute the longitude in [-pi, pi]
    let lon = self.y.atan2(self.x);
    // Conforms to convention: Longitude in [0, 2*PI]
    LonLat::new(
      if lon < 0.0 { 2.0 * PI + lon } else { lon },
      lat
    )
  }

  /// Compute the dot product of this vector with the given vector 
  pub fn scalar(&self, rhs: &XYZ) -> f64 {
    self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
  }
}


/// X, Y coordinates in an image
pub struct ImgXY {
  x: f64,
  y: f64,  
}
impl ImgXY {
  /// New image coordinates fo (x, y).
  pub fn new(x: f64, y: f64) -> Self {
    Self { x, y }
  }

  /// Get the x coordinate
  pub fn x(&self) -> f64 {
    self.x
  }

  /// Get the y coordinate
  pub fn y(&self) -> f64 {
    self.y
  }
}

/*
/// Intermediate World coordinates.
pub struct InterXY {
  x: f64,
  y: f64,
}*/

/// X, Y coordinates in the 2D projection plane
#[derive(Debug, Clone, PartialEq)]
pub struct ProjXY {
  x: f64,
  y: f64,
}
impl ProjXY {
  pub fn new(x: f64, y: f64) -> Self {
    Self { x, y } 
  }

  /// Get the x coordinate
  pub fn x(&self) -> f64 {
    self.x
  }

  /// Get the y coordinate
  pub fn y(&self) -> f64 {
    self.y
  }
}

/// Generic projection (i.e. not necessarily  centered around the vernal point).
pub trait Projection {

  /// Get the projection short name.
  fn short_name(&self) -> &'static str;

  /// Returns the (X, Y) bounds of the projection in the projection plane (`None` if unbounded).
  /// # Remark
  /// * This is not `const` because the bound may depends on projection parameters 
  ///   (see e.g. `CYP` or `CEA`).
  fn bounds(&self) -> &ProjBounds;

  /// Returns `true` if the given point is in the valid projection area.
  fn is_in_valid_proj_area(&self, pos: &ProjXY) -> bool;
  
  /// Project (if possible) from the unit sphere to a projection 2D plane.
  fn proj_xyz(&self, xyz: &XYZ) -> Option<ProjXY>;

  /// Deproject (if possible) from a projection 2D plane to the unit sphere.
  fn unproj_xyz(&self, pos: &ProjXY) -> Option<XYZ>;

  /// Project (if possible) from equatorial coordinates to a 2D projection plane. 
  fn proj_lonlat(&self, lonlat: &LonLat) -> Option<ProjXY> {
    self.proj_xyz(&lonlat.to_xyz())
  }

  /// Deproject (if possible) from a 2D projection plane to equatorial coordinates.
  fn unproj_lonlat(&self, pos: &ProjXY) -> Option<LonLat> {
    self.unproj_xyz(pos).map(|xyz| xyz.to_lonlat())
  }
}

/// The X and Y ranges bounds of a projection in the Euclidean projection plane.
#[derive(Debug, Clone)]
pub struct ProjBounds {
  x: Option<RangeInclusive<f64>>,
  y: Option<RangeInclusive<f64>>,
}

impl ProjBounds {
  
  /// Create a new projection bounds.
  const fn new(x: Option<RangeInclusive<f64>>, y: Option<RangeInclusive<f64>>) -> Self {
    Self { x, y }
  }
  
  /// Returns the bounds of a projection along the x-axis (`None` mean unbounded).
  pub fn x_bounds(&self) -> &Option<RangeInclusive<f64>> {
    &self.x
  }

  /// Returns the bounds of a projection along the y-axis (`None` mean unbounded).
  pub fn y_bounds(&self) -> &Option<RangeInclusive<f64>> {
    &self.y
  }

}


// https://www.aanda.org/articles/aa/full/2002/45/aah3860/aah3860.html
/// Projection centered around the vernal point.
pub trait CanonicalProjection {
  /// Full projection name 
  const NAME: &'static str;
  /// WCS projection name (3 characters)
  const WCS_NAME: &'static str;
  
  /// Returns the (X, Y) bounds of the projection in the projection plane (`None` if unbounded).
  /// # Remark
  /// * This is not `const` because the bound may depends on projection parameters 
  ///   (see e.g. `CYP` or `CEA`).
  fn bounds(&self) -> &ProjBounds;

  /// Returns `true` if the given point is in the valid projection area.
  fn is_in_valid_proj_area(&self, pos: &ProjXY) -> bool {
    // TODO: add the method in each proj to useless computations when the point is valid
    self.unproj(pos).is_some()
  }
  
  /// Project (if possible) from the unit sphere to the canonical projection 2D plane.
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY>;

  /// Deproject (if possible) from the canonical projection 2D plane to the unit sphere.
  fn unproj(&self, pos: &ProjXY) -> Option<XYZ>;
  
}

impl<T: CanonicalProjection> Projection for T {

  fn short_name(&self) -> &'static str {
    Self::WCS_NAME
  }

  fn bounds(&self) -> &ProjBounds { 
    CanonicalProjection::bounds(self)
  }

  fn is_in_valid_proj_area(&self, pos: &ProjXY) -> bool {
    CanonicalProjection::is_in_valid_proj_area(self, pos)
  }

  fn proj_xyz(&self, xyz: &XYZ) -> Option<ProjXY> {
    self.proj(xyz)
  }
  fn unproj_xyz(&self, pos: &ProjXY) -> Option<XYZ> {
    self.unproj(pos)
  }
}

/// Structure performing a rotation (due to non-vernal projection origin)
/// before projecting/after deprojecting.
pub struct CenteredProjection<T: CanonicalProjection> {
  // Parameters of the rotation matrix
  r11: f64, r12: f64, r13: f64,
  r21: f64, r22: f64, r23: f64,
  r31: f64, r32: f64, r33: f64,
  /// Internal projection
  proj: T,
}

impl<T: CanonicalProjection> CenteredProjection<T> {
  
  /// Create a new projection, centered on the vernal point.
  pub fn new(proj: T) -> Self {
    Self {
      r11: 1.0, r12: 0.0, r13: 0.0,
      r21: 0.0, r22: 1.0, r23: 0.0,
      r31: 0.0, r32: 0.0, r33: 1.0,
      proj
    }
  }
  
  pub fn inner_proj(&self) -> &T {
    &self.proj
  }
  
  /// Change the projection center.
  /// # Param
  /// * `lonlat`: new projection center
  pub fn set_proj_center_from_lonlat(&mut self, lonlat: &LonLat) {
    let (sinl, cosl) = lonlat.lon.sin_cos();
    let (sinb, cosb) = lonlat.lat.sin_cos();
    self.r11 = cosl * cosb;  self.r12 = sinl * cosb;  self.r13 = sinb;
    self.r21 = -sinl;        self.r22 =  cosl;        self.r23 =  0.0;
    self.r31 = -cosl * sinb; self.r32 = -sinl * sinb; self.r33 = cosb;
  }

  /// Change the projection center.
  /// # Param
  /// * `xyz`: new projection center
  pub fn set_proj_center_from_xyz(&mut self, xyz: &XYZ) {
    // x = cos(l) * cos(b)
    // y = sin(l) * cos(b)
    // z = sin(b)
    let sinb = xyz.z;
    let cosb = (xyz.x.pow2() + xyz.y.pow2()).sqrt(); // cos = sqrt(1.0 - sin^2)
    let (sinl, cosl) = if cosb == 0.0 {
      (0.0, 1.0)
    } else {
      (xyz.y / cosb,  xyz.x / cosb)
    }; 
    self.r11 = xyz.x;        self.r12 = xyz.y;        self.r13 = xyz.z;
    self.r21 = -sinl;        self.r22 =  cosl;        self.r23 = 0.0;
    self.r31 = -cosl * sinb; self.r32 = -sinl * sinb; self.r33 = cosb;
  }

  /// This describes the computation of a rotation matrix from 3 euler angles
  ///
  /// First rotation is done by a rotation of -lon around z axis (longitude increasing towards the east by default)
  /// Second rotation is done by a rotation of +lat around the y' axis
  /// Third rotation is done by a rotation of +gamma around the x'' axis
  pub fn set_proj_center_from_lonlat_and_positional_angle(&mut self, lonlat: &LonLat, gamma: f64) {
    let (sinl, cosl) = lonlat.lon.sin_cos();
    let (sinb, cosb) = lonlat.lat.sin_cos();
    let (sing, cosg) = gamma.sin_cos();

    self.r11 = cosl * cosb;                       self.r12 = sinl * cosb;                      self.r13 = sinb;
    self.r21 = -sinl * cosg + cosl * sinb * sing; self.r22 = cosl * cosg + sinl * sinb * sing; self.r23 = -cosb * sing;
    self.r31 = -sinl * sing - cosl * sinb * cosg; self.r32 = cosl * sing - sinl * sinb * cosg; self.r33 = cosb * cosg;
  }
}

impl<T: CanonicalProjection> Projection for CenteredProjection<T> {

  fn short_name(&self) -> &'static str {
    self.inner_proj().short_name()
  }

  fn bounds(&self) -> &ProjBounds {
    self.inner_proj().bounds()
  }

  fn is_in_valid_proj_area(&self, pos: &ProjXY) -> bool {
    self.inner_proj().is_in_valid_proj_area(pos)
  }
  
  fn proj_xyz(&self, xyz: &XYZ) -> Option<ProjXY> {
    let rotated_xyz = XYZ::new(
      self.r11 * xyz.x + self.r12 * xyz.y + self.r13 * xyz.z,
      self.r21 * xyz.x + self.r22 * xyz.y + self.r23 * xyz.z,
      self.r31 * xyz.x + self.r32 * xyz.y + self.r33 * xyz.z,
    );
    self.proj.proj(&rotated_xyz) 
  }

  fn unproj_xyz(&self, pos: &ProjXY) -> Option<XYZ> {
    self.proj.unproj(pos)
      .map(|XYZ{x, y, z}|
        XYZ::new(
           self.r11 * x + self.r21 * y + self.r31 * z,
           self.r12 * x + self.r22 * y + self.r32 * z,
           self.r13 * x + self.r23 * y + self.r33 * z
        )
      )
  }
}
