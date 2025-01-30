//! Module containing all implemented pseudocylindrical projections
//! (each projection is in its own sub-module).
//!
//! > Pseudocylindricals projections are like cylindrical projection except that
//! > the parallels of latitude are projected at diminishing lengths towards the
//! > polar regions in order to reduce lateral distorsion there.
//! > (Calabretta and Greisen).

pub mod ait;
pub mod mol;
pub mod par;
pub mod sfl;
