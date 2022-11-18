//! Module containing the structure to convert back on forth 
//! from Image coordinates to Celestial coordinates.

use crate::img2proj::{ImgXY2ProjXY, ProjXY2ImgXY};
use crate::{ImgXY, LonLat, Projection};

/// Structure to convert back on forth from Image coordinates to Celestial coordinates.
pub struct Img2LonLat<P: Projection> {
  img2proj: ImgXY2ProjXY,
  proj2img: ProjXY2ImgXY,
  proj: P,
}

impl<P: Projection> Img2LonLat<P> {
  
  pub fn new(img2proj: ImgXY2ProjXY, proj: P) -> Self {
    let proj2img = img2proj.inverse();
    Self { img2proj, proj2img, proj }
  }
  
  pub fn lonlat2img(&self, lonlat: &LonLat) -> Option<ImgXY> {
    self.proj.proj_lonlat(lonlat).and_then(|xy| self.proj2img.proj2img(&xy))
  }
  
  pub fn img2lonlat(&self, img_pos: &ImgXY) -> Option<LonLat> {
    self.proj.unproj_lonlat(&self.img2proj.img2proj(img_pos))
  }
  
}
