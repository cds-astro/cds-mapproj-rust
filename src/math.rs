//! Utility module defining a few mathematical functions like sinc (sinus cardinal = sin(x) / x).

use std::f64::consts::PI;

pub const HALF_PI: f64 = 0.5 * PI;

pub trait CustomFloat {
  /// `2 * x`
  fn twice(self) -> f64;
  /// `x / 2`
  fn half(self) -> f64;
  /// Square: `x^2`
  fn pow2(self) -> Self;
  /// Cardinal sine, i.e. `sin(x) / x`
  fn sinc(self) -> Self;
  /// Cardinal arcsine, i.e. `asin(x) / x`.  
  fn asinc(self) -> Self;
  /// Hyperbolic inverse tangent, i.e. `tanh-1(x)`.
  /// # Warning
  /// * return `Nan` if the value is in `]-1, 1[`.
  fn atanh(self) -> f64;

}

impl CustomFloat for f64 {
  
  fn twice(self) -> f64 {
    2.0 * self // self + self (I hope the compiler know the simple shift bit to be used for x2)
  }

  fn half(self) -> f64 {
    0.5 * self
  }
  
  fn pow2(self) -> Self {
    self * self
  }
  
  fn sinc(self) -> Self {
    if self.abs() > 1.0e-4 {
      self.sin() /  self
    } else {
      // If a is mall, use Taylor expansion of asin(a) / a
      // a = 1e-4 => a^4 = 1.e-16
      let x2 = self.pow2();
      1.0 - x2 * (1.0 - x2 / 20.0) / 6.0
    }
  }
  
  fn asinc(self) -> Self {
    if self.abs() > 1.0e-4 {
      self.asin() / self
    } else {
      // If a is mall, use Taylor expansion of asin(a) / a
      // a = 1e-4 => a^4 = 1.e-16
      let x2 = self.pow2();
      1.0 + x2 * (1.0 + x2 * 9.0 / 20.0) / 6.0
    }
  }
  
  fn atanh(self) -> f64 {
    0.5 * ((1.0 + self) / (1.0 - self)).ln()
  }
  
}
