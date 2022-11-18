# `mapproj`

Implementation of (a part of) map projections defined in the [FITS World Coordinate System (WCS)](https://fits.gsfc.nasa.gov/fits_wcs.html).


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
For zeithal projections, we thus project on the yz-plane 
using the euclidean coordinates instead of the equatorial coordinates.
Changing the projection center correspond to a simple 3x3 matrix multiplication
of the euclidean coordinates.

This work has been first done internally in 2017 (still by F.-X. Pineau), 
in Java, to support more projections in Aladin (not released yet, not public yet).
F. Ochsenbein previously implemented TAN, STG, SIN, ZEA, ARC, AIT, SFL, MER and CEA
(for lambda=1, i.e. Lambert's projection) in its AstroCoo library.


To Do list
----------

* [X] Add conic projections (`COD`, `COE`, `COO`, `COP`)
* [X] Add cylindrical projections (`CAR`, `CEA`, `CYP`, `MER`)
* [X] Add hybrid projection (`HPX`)
* [X] Add pseudo cylindrical projections (`AIT`, `MOL`, `PAR`, `SFL`)
* [X] Add zenithal projections (`AIR`, `ARC`, `AZP`, `FEYE`, `NCP`, `SIN`, `STG`, `SZP`, `TAN`, `ZEA`, `ZPN`)
* [ ] Add polyconic and pseudoconic projections (`BON, PCO`)?
* [ ] Add quad cube projections (`TSC`, `CSC`, `QSC`)?
* [ ] Check and possibly add constants to match WCS projection bounds
* [X] Support `CRPIX` + `CD` convention
* [X] Support `CRPIX` + `PC` + `CDELT` convention
* [X] Support `CRPIX` + `CROTA` + `CDELT` convention
* [ ] Add support for LONPOLE?
* [ ] Test and complete SIP
* [ ] Check and add to git the pdf document containing computational details
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


