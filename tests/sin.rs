use mapproj::{
  CenteredProjection, ImgXY, LonLat,
  img2celestial::Img2Celestial,
  img2proj::WcsImgXY2ProjXY,
  zenithal::sin::Sin,
};

#[test]
fn test_sin() {
  // Test fomr the foloowing FITS cards
  /*
  CTYPE1  = 'RA---SIN' 
  CTYPE2  = 'DEC--SIN' 
  CRPIX1  =      382.00001513958 / reference pixel in X coordinate
  CRPIX2  =     389.500015437603 / reference pixel in Y coordinate
  CRVAL1  =        183.914583333 / RA of reference position (degrees)
  CRVAL2  =              36.3275 / DEC of reference position (degrees)
  WCSDIM  =                    2
  CD1_1   =  -2.7777777349544E-4
  CD2_2   =  2.77777773495436E-4
  CDELT1  =  -2.7777777349544E-4 // Redundancy with CD, we ignore this
  CDELT2  =  2.77777773495436E-4 // Redundancy with CD, we ignore this
  */
  
  // Define constants
  let crpix1 =  382.00001513958_f64;
  let crpix2 = 389.500015437603_f64;
  
  let crval1 = 183.914583333_f64;
  let crval2 = 36.3275_f64;
  
  let cd11 = -2.7777777349544e-4_f64;
  let cd22 =  2.77777773495436e-4_f64;
  
  // Set the projection
  let mut proj = CenteredProjection::new(Sin::default());
  let proj_center = LonLat::new(crval1.to_radians(), crval2.to_radians());
  proj.set_proj_center_from_lonlat(&proj_center);
  let img2proj = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, 0.0, 0.0, cd22);
  let img2lonlat = Img2Celestial::new(img2proj, proj);
  // We could have set the projection center here instead of previously:
  //   img2lonlat.set_proj_center_from_lonlat(proj_center);
  
  // Use to project, unproject coordinates:
  // - we choose on purpose position in the image of the projection center
  let img_coo_input = ImgXY::new(382.00001513958, 389.500015437603);
  let lonlat = img2lonlat.img2lonlat(&img_coo_input).unwrap();
  assert!((lonlat.lon() - proj_center.lon()).abs() < 1e-14);
  assert!((lonlat.lat() - proj_center.lat()).abs() < 1e-14);
  let img_coo_input = img2lonlat.lonlat2img(&lonlat).unwrap();
  assert!((img_coo_input.x() - img_coo_input.x()).abs() < 1e-14);
  assert!((img_coo_input.y() - img_coo_input.y()).abs() < 1e-14);
  
}