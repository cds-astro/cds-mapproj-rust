//! Test to verify SIP coefficient polynomial evaluation

use mapproj::sip::SipCoeff;

#[test]
fn test_sip_poly_evaluation() {
  // Create a simple SIP polynomial for degree 2
  // Coefficients in standard SIP ordering: (0,0), (1,0), (2,0), (0,1), (1,1), (0,2)
  // f(u,v) = c00 + c10*u + c20*u^2 + c01*v + c11*u*v + c02*v^2

  let coeffs = vec![
    1.0, // (0,0)
    2.0, // (1,0)
    3.0, // (2,0)
    4.0, // (0,1)
    5.0, // (1,1)
    6.0, // (0,2)
  ];

  let sip = SipCoeff::new(coeffs.into_boxed_slice());

  // Test at u=2.0, v=3.0
  // Expected: 1 + 2*2 + 3*4 + 4*3 + 5*2*3 + 6*9
  //         = 1 + 4 + 12 + 12 + 30 + 54
  //         = 113
  let result = sip.p(2.0, 3.0);
  let expected = 113.0;

  assert!(
    (result - expected).abs() < 1e-10,
    "Expected {expected}, got {result}"
  );

  // Test at u=0, v=0 (should be just c00)
  let result_origin = sip.p(0.0, 0.0);
  assert!((result_origin - 1.0).abs() < 1e-10);
}

#[test]
fn test_sip_poly_degree4() {
  // Test with actual A coefficients from 74721b067-w2-int-1b.fits
  // Order in standard SIP: for degree 4, we have 15 coefficients
  // (0,0), (1,0), (2,0), (3,0), (4,0), (0,1), (1,1), (2,1), (3,1), (0,2), (1,2), (2,2), (0,3), (1,3), (0,4)

  let a_coeffs = vec![
    0.0,             // A_0_0 (always 0 in SIP)
    -6.45286001e-4,  // A_1_0
    -1.27125200e-6,  // A_2_0
    3.10842100e-9,   // A_3_0
    1.49894000e-12,  // A_4_0
    -1.20498000e-4,  // A_0_1
    9.47114200e-7,   // A_1_1
    -1.16878700e-9,  // A_2_1
    -5.99954500e-13, // A_3_1
    -8.08115600e-6,  // A_0_2
    4.64597500e-10,  // A_1_2
    1.50566800e-12,  // A_2_2
    -3.12543100e-10, // A_0_3
    5.95993800e-14,  // A_1_3
    1.33778800e-13,  // A_0_4
  ];

  let sip = SipCoeff::new(a_coeffs.into_boxed_slice());

  // Test at u=50.63, v=223.08 (from our test pixel)
  let u = 50.63;
  let v = 223.08;
  let result = sip.p(u, v);

  // Expected from Python calculation: -0.456283906403315
  let expected = -0.456283906403315;

  assert!(
    (result - expected).abs() < 1e-6,
    "Expected {expected:.15}, got {result:.15}, error = {:.3e}",
    (result - expected).abs()
  );
}
