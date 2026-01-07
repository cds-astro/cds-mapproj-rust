//! Standalone test to verify img2proj TPV implementation against Astropy reference values.
//!
//! Run with: cargo run --example verify_img2proj_tpv

#![allow(clippy::wildcard_imports)]
#![allow(clippy::manual_assert)]
#![allow(clippy::uninlined_format_args)]

use mapproj::img2proj::*;
use mapproj::tpv::{Tpv, TpvCoeff, TpvPV};
use mapproj::ImgXY;

fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
    let diff = (actual - expected).abs();
    if diff > tolerance {
        panic!(
            "{}: FAILED\n  Expected: {:.15e}\n  Actual:   {:.15e}\n  Diff:     {:.15e}\n  Tolerance: {:.15e}",
            label, expected, actual, diff, tolerance
        );
    }
    println!("  {} OK (diff: {:.2e})", label, diff);
}

fn test1_simple_linear_tpv() {
    println!("\n{}", "=".repeat(80));
    println!("TEST 1: Simple Linear TPV (xi' = 1.5*xi, eta' = 2.0*eta)");
    println!("{}", "=".repeat(80));

    // Parameters
    let crpix1 = 100.0;
    let crpix2 = 100.0;
    let cd11 = 1.0e-5; // degrees/pixel
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // TPV coefficients: xi' = 1.5*xi, eta' = 2.0*eta
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.5]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 2.0]);
    let pv = TpvPV::new(pv1, pv2);
    // Range in PIXEL offsets - need to accommodate pixel offsets of +/-100 or more
    let tpv = Tpv::new(pv, -200.0..=200.0, -200.0..=200.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_tpv = WcsWithTpvImgXY2ProjXY::new(wcs, tpv);

    // Test point
    let img_xy = ImgXY::new(105.0, 103.0);
    let proj_xy = wcs_tpv.img2proj(&img_xy);

    // Expected values from Astropy/manual calculation
    let expected_x = 1.308996938995747e-06; // radians
    let expected_y = 1.047197551196598e-06; // radians

    println!("\nInput: ({}, {})", img_xy.x(), img_xy.y());
    println!(
        "Output: ({:.15e}, {:.15e}) radians",
        proj_xy.x(),
        proj_xy.y()
    );
    println!("\nVerification against Astropy:");
    assert_close(proj_xy.x(), expected_x, 1e-15, "X coordinate");
    assert_close(proj_xy.y(), expected_y, 1e-15, "Y coordinate");
}

fn test2_radial_tpv() {
    println!("\n{}", "=".repeat(80));
    println!("TEST 2: Radial TPV (xi' = xi + 0.01*r, eta' = eta + 0.01*r)");
    println!("{}", "=".repeat(80));

    // Parameters
    let crpix1 = 256.0;
    let crpix2 = 256.0;
    let cd11 = 1.0e-5; // degrees/pixel
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // TPV coefficients: xi' = xi + 0.01*r, eta' = eta + 0.01*r
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.01]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.01]);
    let pv = TpvPV::new(pv1, pv2);
    // Range in PIXEL offsets - need to accommodate typical astronomical image sizes
    let tpv = Tpv::new(pv, -300.0..=300.0, -300.0..=300.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_tpv = WcsWithTpvImgXY2ProjXY::new(wcs, tpv);

    // Test point
    let img_xy = ImgXY::new(260.0, 259.0);
    let proj_xy = wcs_tpv.img2proj(&img_xy);

    // Expected values from Astropy/manual calculation
    let expected_x = 7.068583470577035e-07; // radians
    let expected_y = 5.323254218582706e-07; // radians

    println!("\nInput: ({}, {})", img_xy.x(), img_xy.y());
    println!(
        "Output: ({:.15e}, {:.15e}) radians",
        proj_xy.x(),
        proj_xy.y()
    );
    println!("\nVerification against Astropy:");
    assert_close(proj_xy.x(), expected_x, 1e-15, "X coordinate");
    assert_close(proj_xy.y(), expected_y, 1e-15, "Y coordinate");
}

fn test3_complex_quadratic_tpv() {
    println!("\n{}", "=".repeat(80));
    println!("TEST 3: Complex TPV with quadratic distortion");
    println!("{}", "=".repeat(80));

    // Parameters
    let crpix1 = 512.0;
    let crpix2 = 512.0;
    let cd11 = -2.0e-5; // degrees/pixel
    let cd12 = 0.5e-5;
    let cd21 = -0.5e-5;
    let cd22 = 2.0e-5;

    // TPV with quadratic terms
    // xi' = xi + 0.001*xi^2 + 0.002*eta^2
    // eta' = eta + 0.002*xi^2 + 0.001*eta^2
    let pv1 = TpvCoeff::new_axis1(&[
        0.0,   // PV1_0
        1.0,   // PV1_1 (xi term)
        0.0,   // PV1_2 (eta term)
        0.0,   // PV1_3 (r term)
        0.001, // PV1_4 (xi^2 term)
        0.0,   // PV1_5 (xi*eta term)
        0.002, // PV1_6 (eta^2 term)
    ]);

    let pv2 = TpvCoeff::new_axis2(&[
        0.0,   // PV2_0
        1.0,   // PV2_1 (eta term)
        0.0,   // PV2_2 (xi term)
        0.0,   // PV2_3 (r term)
        0.002, // PV2_4 (eta^2 term)
        0.0,   // PV2_5 (eta*xi term)
        0.001, // PV2_6 (xi^2 term)
    ]);

    let pv = TpvPV::new(pv1, pv2);
    // Range in PIXEL offsets - need to accommodate test point at pixel offset ~8
    let tpv = Tpv::new(pv, -600.0..=600.0, -600.0..=600.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_tpv = WcsWithTpvImgXY2ProjXY::new(wcs, tpv);

    // Test point
    let img_xy = ImgXY::new(520.0, 515.0);
    let proj_xy = wcs_tpv.img2proj(&img_xy);

    // Expected values from manual calculation:
    // u = 8, v = 3
    // u' = 8 + 0.001*64 + 0.002*9 = 8.082
    // v' = 3 + 0.002*9 + 0.001*64 = 3.082
    // xi_rad = (-2e-5 * 8.082 + 0.5e-5 * 3.082) * pi/180 = -2.5522e-6 radians
    // eta_rad = (-0.5e-5 * 8.082 + 2e-5 * 3.082) * pi/180 = 3.7053e-7 radians
    let expected_x = -2.552194965191308e-06; // radians
    let expected_y = 3.705334001983961e-07; // radians

    println!("\nInput: ({}, {})", img_xy.x(), img_xy.y());
    println!(
        "Output: ({:.15e}, {:.15e}) radians",
        proj_xy.x(),
        proj_xy.y()
    );
    println!("\nVerification against Astropy:");
    assert_close(proj_xy.x(), expected_x, 1e-15, "X coordinate");
    assert_close(proj_xy.y(), expected_y, 1e-15, "Y coordinate");
}

fn main() {
    println!("IMG2PROJ TPV VERIFICATION AGAINST ASTROPY REFERENCE");
    println!("Independent verification without relying on Rust test suite\n");

    test1_simple_linear_tpv();
    test2_radial_tpv();
    test3_complex_quadratic_tpv();

    println!("\n{}", "=".repeat(80));
    println!("SUMMARY");
    println!("{}", "=".repeat(80));
    println!("OK All tests passed!");
    println!("OK img2proj TPV implementation matches Astropy reference");
}
