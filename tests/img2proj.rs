use mapproj::{
  ImgXY,
  img2proj::{ImgXY2ProjXY, ProjXY2ImgXY, WcsImgXY2ProjXY},
};

fn assert_round_trip(img2proj: &WcsImgXY2ProjXY) {
  let proj2img = img2proj.inverse();
  for (x, y) in [(8.25, 14.5), (0.0, 0.0), (-3.5, 21.0)] {
    let xy = img2proj.img2proj(&ImgXY::new(x, y));
    let img = proj2img.proj2img(&xy).unwrap();
    assert!((img.x() - x).abs() < 1e-12, "x: {} != {}", img.x(), x);
    assert!((img.y() - y).abs() < 1e-12, "y: {} != {}", img.y(), y);
  }
}

#[test]
fn test_wcs_inverse_non_symmetric_cd() {
  // CD1_1 = 2, CD1_2 = 1, CD2_1 = 0, CD2_2 = 3
  assert_round_trip(&WcsImgXY2ProjXY::from_cd(7.25, 12.5, 2.0, 1.0, 0.0, 3.0));
  assert_round_trip(&WcsImgXY2ProjXY::from_cd(7.25, 12.5, -0.07, -0.0189, -0.0078, 0.057));
  assert_round_trip(&WcsImgXY2ProjXY::from_pc(7.25, 12.5, 1.0, 0.27, -0.13, 0.95, 2.0, 3.0));
  assert_round_trip(&WcsImgXY2ProjXY::from_cr(7.25, 12.5, 30.0, 2.0, 3.0));
}

#[test]
fn test_wcs_crota2() {
  // CDELTi + CROTA2 as converted by wcslib: CD1_2 = -CDELT2 sin(CROTA2), CD2_1 = CDELT1 sin(CROTA2)
  let (crota2, cdelt1, cdelt2) = (30.0_f64, 2.0, 3.0);
  let (s, c) = crota2.to_radians().sin_cos();
  let from_cr = WcsImgXY2ProjXY::from_cr(0.0, 0.0, crota2, cdelt1, cdelt2);
  let from_cd = WcsImgXY2ProjXY::from_cd(0.0, 0.0, cdelt1 * c, -cdelt2 * s, cdelt1 * s, cdelt2 * c);
  for (x, y) in [(1.0, 0.0), (0.0, 1.0), (-2.5, 4.0)] {
    let a = from_cr.img2proj(&ImgXY::new(x, y));
    let b = from_cd.img2proj(&ImgXY::new(x, y));
    assert!((a.x() - b.x()).abs() < 1e-14, "x: {} != {}", a.x(), b.x());
    assert!((a.y() - b.y()).abs() < 1e-14, "y: {} != {}", a.y(), b.y());
  }
}
