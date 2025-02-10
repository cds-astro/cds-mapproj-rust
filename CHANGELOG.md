# `mapproj` Change Log


## 0.4.0

Released 2025-02-10

### Add

* Add a constructor for defining a rotation matrix from a lonlat + position angle


## 0.3.0

Released 2022-12-13

### Fix

* Conic COE proj (x instead of z in projection) 
* Conic unproj arctan (`r` parameter was important because of possible negative sign)

### Changed

* Add `bounds` and  `is_in_valid_proj_area` methods projections
* In `unproj` check if the (X, Y) cooridnates are in the projection plane valid area
* `ProjXY2ImgXY` now an trait with implementations.
* Support SIP in a different struct implementing `ProjXY2ImgXY`
* Add `BasicImgXY2ProjXY` and `ReversedEastPngImgXY2ProjXY`
* Conic proj default parameters

--------------------------------------------------------------------------------


## 0.2.0

Released 2022-11-18

### Changed

* rename `Img2LonLat` in `Img2Celestial` and add the `XYZ` in input/output
  in addition to `LonLat`

### Add

* Add `sin` example in test and in README
* Add public method `x()`, `y()`, `z()`
* Add `set_proj_center_from_xx` in `Img2Celestial`

--------------------------------------------------------------------------------


## 0.1.0

Released 2022-11-18


