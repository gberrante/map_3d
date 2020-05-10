# Rust Map3d
![Rust](https://github.com/gberrante/map_3d/workflows/Rust/badge.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/cab9b019a46644f59dce4b3b21d5404a)](https://www.codacy.com/manual/errante.gianni/map_3d?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=gberrante/map_3d&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/gberrante/map_3d/branch/master/graph/badge.svg)](https://codecov.io/gh/gberrante/map_3d)

This is a Rust library for geographic coordinate frame conversion. The implementation is similar to  [Pymap3d](https://github.com/geospace-code/pymap3d). All the functions are implemented in `f64` precision. 

No external libraries are needed.

The default units are:

- Radians [rad] for angular variables
- Meters  [m] for linear variables
- [Greenwich Sidereal Time](https://www.cfa.harvard.edu/~jzhao/times.html)  [GST] for date and time

The default reference ellipsoid is the [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System#A_new_World_Geodetic_System:_WGS_84)

List of coordinates systems implemented in the functions:

- [Geodetic coordinate system](https://en.wikipedia.org/wiki/Geographic_coordinate_system) (GEODETIC)

- [Earth-Centered Earth-fixed](https://en.wikipedia.org/wiki/ECEF) (ECEF)
- [Earth-Centered Inertial](https://en.wikipedia.org/wiki/Earth-centered_inertial) (ECI) 
- [Local spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system#In_geography) (AER)
- [Local tangent plane coordinate system - East-North-Up](https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates) (ENU)



## To-Do List:

- implement function for UTC conversion to GST
- implement functions for geodetic distances calculation
- implement functions for North-East-Down coordinate system
- implement haversine functions
- implement functions for right ascension and declination conversions

- create HTML documentation