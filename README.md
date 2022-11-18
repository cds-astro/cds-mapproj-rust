# `mapproj`

Implementation of (a part of) map projections defined in the [FITS World Coordinate System (WCS)](https://fits.gsfc.nasa.gov/fits_wcs.html).

[![](https://img.shields.io/crates/v/mapproj.svg)](https://crates.io/crates/mapproj)
[![](https://img.shields.io/crates/d/mapproj.svg)](https://crates.io/crates/mapporj)
[![API Documentation on docs.rs](https://docs.rs/mapproj/badge.svg)](https://docs.rs/mapproj/)


Purpose
-------

Mainly to add projections and to support the display of FITS images in Aladin Lite V3.


Warning
-------

This library:
* does not support FITS reading/parsing and keywords analysis;
* is still in an early stage (see To Do List) and will evolve when included in Aladin Lite V3.


Details
-------

Contrary to the [WCS paper](https://arxiv.org/pdf/astro-ph/0207413.pdf), 
and following a previous work by F. Ochsenbein, we use the vernal point
`(1, 0, 0)` as default projection center/origin. 
For zenithal projections, we thus project on the yz-plane 
using the euclidean coordinates instead of the equatorial coordinates.
Changing the projection center corresponds to a simple 3x3 matrix multiplication
of the euclidean coordinates.

This work has been first done internally in 2017 (still by F.-X. Pineau), 
in Java, to support more projections in Aladin Desktop (not released yet, not public yet).
F. Ochsenbein previously implemented TAN, STG, SIN, ZEA, ARC, AIT, SFL, MER and CEA
(for lambda=1, i.e. Lambert's projection) in its AstroCoo library.


Example
-------

Given the following FITS cards:
```bash
CTYPE1  = 'RA---SIN'                                               
CTYPE2  = 'DEC--SIN'                            
CRPIX1  =      382.00001513958 / reference pixel in X coordinate  
CRPIX2  =     389.500015437603 / reference pixel in Y coordinate       
CRVAL1  =        183.914583333 / RA of reference position (degrees)        
CRVAL2  =              36.3275 / DEC of reference position (degrees)      
WCSDIM  =                    2                            
CD1_1   =  -2.7777777349544E-4                                  
CD2_2   =  2.77777773495436E-4                                                                               
CDELT1  =  -2.7777777349544E-4 / Redundancy with CD1_1, we ignore it      
CDELT2  =  2.77777773495436E-4 / Redundancy with CD2_2, we ignore it
```

```rust
// Imports
use mapproj::{
   CenteredProjection, ImgXY, LonLat
   img2lonlat::Img2LonLat,
   img2proj::ImgXY2ProjXY,
   zenithal::sin::Sin,
};

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
let img2proj = ImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, 0.0, 0.0, cd22);
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
```


To Do list
----------

* [X] Add conic projections (`COD`, `COE`, `COO`, `COP`)
* [X] Add cylindrical projections (`CAR`, `CEA`, `CYP`, `MER`)
* [X] Add hybrid projection (`HPX`)
* [X] Add pseudo cylindrical projections (`AIT`, `MOL`, `PAR`, `SFL`)
* [X] Add zenithal projections (`AIR`, `ARC`, `AZP`, `FEYE`, `NCP`, `SIN`, `STG`, `SZP`, `TAN`, `ZEA`, `ZPN`)
* [ ] Add polyconic and pseudoconic projections (`BON, PCO`)?
* [ ] Add quad cube projections (`TSC`, `CSC`, `QSC`)?
* [ ] Check and possibly document constants to be added to match WCS projection bounds
* [X] Support `CRPIX` + `CD` convention
* [X] Support `CRPIX` + `PC` + `CDELT` convention
* [X] Support `CRPIX` + `CROTA` + `CDELT` convention
* [ ] Add support for LONPOLE?
* [ ] Test and complete SIP
* [X] Add to git the pdf document containing computational details
* [ ] Check, fix typo, enrich the pdf document containing computational details
* [ ] Add generation of projection files and plots (like in the Java lib)


License
-------

Like most projects in Rust, this project is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or
   http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or
   http://opensource.org/licenses/MIT)

at your option.


Contribution
------------

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this project by you, as defined in the Apache-2.0 license,
shall be dual licensed as above, without any additional terms or conditions.


