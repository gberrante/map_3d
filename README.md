# Rust Map3d
[![crates.io](https://img.shields.io/crates/v/map_3d.svg)](https://crates.io/crates/map_3d)
![Rust](https://github.com/gberrante/map_3d/workflows/Rust/badge.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![codecov](https://codecov.io/gh/gberrante/map_3d/branch/master/graph/badge.svg)](https://codecov.io/gh/gberrante/map_3d)

This is a Rust library for geographic coordinate frame conversion. The implementation is similar to  [Pymap3d](https://github.com/geospace-code/pymap3d). All the functions are implemented in `f64` precision. 

Live demo: [map 3d live demo](https://rustmap-3d.firebaseapp.com/)

No external dependencies

The default units are:

- Radians [rad] for angular variables
- Meters  [m] for linear variables
- [Greenwich Sidereal Time](https://www.cfa.harvard.edu/~jzhao/times.html)  [GST] for date and time

We support 4 reference ellipsoids at the moment, and
[WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System#A_new_World_Geodetic_System:_WGS_84).
is the default one, obtained with `Reference::default()`.

List of coordinates systems implemented in the functions:

- [Geodetic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system) (GEODETIC)
- [Earth-Centered Earth-fixed](https://en.wikipedia.org/wiki/ECEF) (ECEF)
- [Earth-Centered Inertial](https://en.wikipedia.org/wiki/Earth-centered_inertial) (ECI) 
- [Local spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system#In_geography) (AER)
- [Local tangent plane coordinate system - East-North-Up](https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates) (ENU)
- [Local tangent plane coordinate system - North-East-Down](https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates) (NED)

Additional functions:

- Radians to Degrees and Degrees to Radians
- UTC time conversion to GST
- 3x3 Matrix - 3x1 column multiplication
- 3x3 Matrix transpose
- f64 round towards zero 
- projected distance (`Haversine` formula) between two coordinates (lat, lon, in decimal degrees)

## To-Do List:

- implement functions for right ascension and declination conversions
