use std::f64::consts::PI;

use mapproj::{
  self, CustomFloat,
  LonLat, ProjXY,
  Projection,
  conic::{
    cod::Cod,
    coe::Coe,
    coo::Coo,
    cop::Cop,
  },
  cylindrical::{
    car::Car,
    cea::Cea,
    cyp::Cyp,
    mer::Mer
  },
  hybrid::hpx::Hpx,
  pseudocyl::{
    ait::Ait,
    mol::Mol,
    par::Par,
    sfl::Sfl,
  },
  zenithal::{
    air::Air,
    arc::Arc,
    azp::Azp,
    feye::Feye,
    ncp::Ncp,
    sin::Sin,
    stg::Stg,
    szp::Szp,
    tan::Tan,
    zea::Zea,
    zpn::Zpn,
  }
};

struct Point {
  color: f64,
  coo: LonLat,
  xy: ProjXY,
}

impl Point {

  fn new(lonlat: LonLat) -> Self {
    let color = (90.0 + lonlat.lat().to_degrees()) / 180.0;
    Self {
      color,
      coo: lonlat,
      xy: ProjXY::new(0.0, 0.0)
    }
  }
  
}

fn gen_coo() -> Vec<LonLat> {
  let mut coos: Vec<LonLat> = Vec::with_capacity(1800 * 12 + 3600 * 11);
  // Meridians (lon = cte): 0, 30, 60, 90, ... 360
  for lon in (0..=360).step_by(30) { // x12
    for lat in -900..=900 { // 180 * 10 = 1800 points
      let mut lon_rad = (lon as f64).to_radians();
      if lon_rad >= PI.twice() {
        lon_rad = PI.twice() - 1.0e-14;
      }
      let lat_rad = ((lat as f64) / 10.0).to_radians();
      coos.push(LonLat::new(lon_rad, lat_rad))
    }
  }
  // Circles of latitude (lat = cte): 0, 15, 30, 45, 60, 75 // (+-15 1+2x5) 
  for lat in (-90..=90).step_by(15) {
    for lon in 0..=3600 {
      let mut lon_rad = ((lon as f64) / 10.0).to_radians();
      if lon_rad >= PI.twice() {
        lon_rad = PI.twice() - 1.0e-14;
      }
      let lat_rad = (lat as f64).to_radians();
      coos.push(LonLat::new(lon_rad, lat_rad))
    }
  }
  // 90 - epislon (10^-10 => 0.3 microarcsec)
  let lat_rad = 90.0_f64.to_radians() - 1.0e-10;
  for lon in 0..=3600 {
    let mut lon_rad = ((lon as f64) / 10.0).to_radians();
    if lon_rad >= PI.twice() {
      lon_rad = PI.twice() - 1e-14;
    }
    coos.push(LonLat::new(lon_rad, lat_rad))
  }
  // -90 + epsilon
  let lat_rad = -90.0_f64.to_radians() + 1.0e-10;
  for lon in 0..=3600 {
    let mut lon_rad = ((lon as f64) / 10.0).to_radians();
    if lon_rad >= PI.twice() {
      lon_rad = PI.twice() - 1.0e-14;
    }
    coos.push(LonLat::new(lon_rad, lat_rad))
  }
  coos
}

fn gen_points() -> Vec<Point> {
  gen_coo()
    .drain(..)
    .map(Point::new)
    .collect()
}

/*
fn build_file<T: Projection>(proj: T, filename: &str) {
  let points = gen_points();
  for point in points {
    proj
  }
}*/

fn test_canonical_back_and_forth<T: Projection>(proj: T, precision_rad: f64) {
  let coos = gen_coo();
  let tot = coos.len();
  let mut i = 0;
  let mut j = 0;
  for coo in coos {
    if let Some(xy) = proj.proj_lonlat(&coo) {
      j += 1;
      if let Some(coo2) = proj.unproj_lonlat(&xy) {
        let ang_dist = coo.haversine_dist(&coo2);
        assert!(
          ang_dist < precision_rad, 
          "i: {}; {} arcsec; x: {}; y: {}; lon: {}; lat: {}, lon2: {}, lat2: {}; {:?}", 
          i, ang_dist.to_degrees() * 3600.0, xy.x(), xy.y(),
          coo.lon().to_degrees(), coo.lat().to_degrees(),
          coo2.lon().to_degrees(), coo2.lat().to_degrees(),
          coo.to_xyz()
        );
        i += 1;
      }
    }
  }
  println!("{}: {}-{}/{}", proj.short_name(), j, i, tot);
}

#[test]
fn test_canonicals_back_and_forth() {
  let mas_in_rad = (1.0 / 3600_000.0_f64).to_radians();
  // Conic
  test_canonical_back_and_forth(Cod::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Coe::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Coo::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Cop::new(), mas_in_rad / 1000.0);        //  1 uas
  // Cylindrical
  test_canonical_back_and_forth(Car::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Cea::new(), 30.0 * mas_in_rad / 1000.0); // 30 uas => due to z = sin(lat)
  test_canonical_back_and_forth(Cyp::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Mer::new(), mas_in_rad / 1000.0);        //  1 uas
  // Hybrid
  test_canonical_back_and_forth(Hpx::new(), mas_in_rad / 1000.0);        //  1 uas
  // Pseudo-Cylindrical
  test_canonical_back_and_forth(Ait::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Mol::new(), 30.0 * mas_in_rad / 1000.0); // 30 uas => also due to operations on z = sin(lat)
  test_canonical_back_and_forth(Par::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Sfl::new(), mas_in_rad / 1000.0);        //  1 uas
  // Zenithal
  test_canonical_back_and_forth(Air::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Arc::new(), 5.0 * mas_in_rad);           //  5 mas => bad precision for lon=90, lat near -90
  test_canonical_back_and_forth(Azp::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Feye::new(), 5.0 * mas_in_rad);          //  5 mas => same as Arc
  test_canonical_back_and_forth(Ncp::new(), 5.0 * mas_in_rad);           //  5 mas => same as Arc
  test_canonical_back_and_forth(Sin::new(), 5.0 * mas_in_rad);           //  5 mas
  test_canonical_back_and_forth(Stg::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Szp::new(), 10000.0 * mas_in_rad);       //  10 arcsec!!
  test_canonical_back_and_forth(Tan::new(), mas_in_rad / 1000.0);        //  1 uas
  test_canonical_back_and_forth(Zea::new(), 2.0 * mas_in_rad);           //  2 mas
  test_canonical_back_and_forth(Zpn::from_params(vec![0.0, 1.0, 0.0, -50.0]).unwrap(), 10.0 * mas_in_rad); // 10 mas
  test_canonical_back_and_forth(Zpn::from_params(vec![0.050, 0.975, -0.807, 0.337, -0.065, 0.010, 0.003, -0.001]).unwrap(), 10.0 * mas_in_rad); // 10 mas
}