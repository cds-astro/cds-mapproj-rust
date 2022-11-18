//! Module containing all implemented conic projections
//! (each projection is in its own sub-module).

pub mod cod;
pub mod coe;
pub mod coo;
pub mod cop;

#[derive(Debug, Clone)]
struct Conic {
  /// Parameters `ThetaA` (in radians)
  ta: f64,
  /// Parameter `Nu` (in radians)
  nu: f64,
  /// Pre-computed quantity `ThetaA - Nu`
  theta1: f64,
  /// Pre-computed quantity `ThetaA + Nu`
  theta2: f64,
  /// Pre-computed quantity telling if `ThetaA` is positive or negative
  negative_ta: bool,
}

impl Conic {
  
  ///
  /// # Params
  /// * `theta_a`: in radians
  /// * `nu`: in radians
  fn from_params(theta_a: f64, nu: f64) -> Self {
    let ta = theta_a;
    let theta1 = ta - nu;
    let theta2 = ta + nu;
    let negative_ta = ta < 0.0;
    Self {
      ta, nu,
      theta1, theta2, negative_ta
    }
  }

}