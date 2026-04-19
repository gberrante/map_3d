//! map_3d: geographic coordinates conversions

/// Returns the tuple (x,y,z) of coordinates in the ECEF system
/// with desired reference frame.
///
/// ## Inputs:
/// - lat = latitude [rad]
/// - lon = longitude [rad]
/// - alt = altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
pub fn geodetic2ecef(lat: f64, lon: f64, alt: f64, r_ellips: Ellipsoid) -> (f64, f64, f64) {
    let n = get_radius_normal(lat, r_ellips);
    let (major, minor, _, _) = r_ellips.parameters();

    let x = (n + alt) * lat.cos() * lon.cos();
    let y = (n + alt) * lat.cos() * lon.sin();
    let z = (n * (minor / major) * (minor / major) + alt) * lat.sin();

    (x, y, z)
}

/// Returns the tuple (azimuth,elevation,slant range) of coordinates in the AER system
///
/// ## Inputs:
/// - lat = input latitude [rad]
/// - lon = input longitude [rad]
/// - alt = input altitude [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - az = azimuth angle [rad] of input geodetic location from reference geodetic location
/// - el = elevation angle [rad] of input geodetic location from reference geodetic location
/// - slant_range = slant range [m] of input geodetic location from reference geodetic location
pub fn geodetic2aer(
    lat: f64,
    lon: f64,
    alt: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (e, n, u) = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, r_ellips);
    let (az, el, slant_range) = enu2aer(e, n, u);

    (az, el, slant_range)
}

/// Returns the tuple (east,north,up) of coordinates in the ENU system
///
///  ## Inputs:
/// - lat = input latitude [rad]
/// - lon = input longitude [rad]
/// - alt = input altitude [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - e = east coordinate [m] of input geodetic location from reference geodetic location
/// - n = north coordinate [m] of input geodetic location from reference geodetic location
/// - u = up coordinate [m] of input geodetic location from reference geodetic location
pub fn geodetic2enu(
    lat: f64,
    lon: f64,
    alt: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x1, y1, z1) = geodetic2ecef(lat, lon, alt, r_ellips);
    let (x2, y2, z2) = geodetic2ecef(lat0, lon0, alt0, r_ellips);

    let (e, n, u) = uvw2enu(x1 - x2, y1 - y2, z1 - z2, lat0, lon0);

    (e, n, u)
}

/// Returns the tuple (north,east,down) of coordinates in the NED system
///
///  ## Inputs:
/// - lat = input latitude [rad]
/// - lon = input longitude [rad]
/// - alt = input altitude [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - n = north coordinate [m] of input geodetic location from reference geodetic location
/// - e = east coordinate [m] of input geodetic location from reference geodetic location
/// - d = down coordinate [m] of input geodetic location from reference geodetic location
pub fn geodetic2ned(
    lat: f64,
    lon: f64,
    alt: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let enu = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, r_ellips);
    (enu.1, enu.0, -enu.2)
}

/// Returns the tuple (x,y,z) of coordinates in the ECEF system
///
/// ## Inputs:
/// - az = azimuth angle [rad]
/// - el = elevation angle [rad]
/// - slant_range = slant range [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
pub fn aer2ecef(
    az: f64,
    el: f64,
    slant_range: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x0, y0, z0) = geodetic2ecef(lat0, lon0, alt0, r_ellips);
    let (e, n, u) = aer2enu(az, el, slant_range);
    let (dx, dy, dz) = enu2uvw(e, n, u, lat0, lon0);
    (x0 + dx, y0 + dy, z0 + dz)
}

/// Returns the tuple (east,north,up) of coordinates in the ENU system
///
/// ## Inputs:
/// - az = azimuth angle [rad]
/// - el = elevation angle [rad]
/// - slant_range = slant range [m]
///
/// ## Outputs:
/// - e = east coordinate [m] of input location from reference geodetic location
/// - n = north coordinate [m] of input location from reference geodetic location
/// - u = up coordinate [m] of input location from reference geodetic location
pub fn aer2enu(az: f64, el: f64, slant_range: f64) -> (f64, f64, f64) {
    let r = slant_range * el.cos();
    (r * az.sin(), r * az.cos(), slant_range * el.sin())
}

/// Returns the tuple (x,y,z) of coordinates in the ECI system
///
/// ## Inputs:
/// - az = azimuth angle [rad]
/// - el = elevation angle [rad]
/// - slant_range = slant range [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - x = x ECI coordinate [m]
/// - y = y ECI coordinate [m]
/// - z = z ECI coordinate [m]
pub fn aer2eci(
    gst: f64,
    az: f64,
    el: f64,
    slant_range: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x1, y1, z1) = aer2ecef(az, el, slant_range, lat0, lon0, alt0, r_ellips);
    ecef2eci(gst, x1, y1, z1)
}

/// Returns the tuple (north,east,down) of coordinates in the NED system
///
/// ## Inputs:
/// - az = azimuth angle [rad]
/// - el = elevation angle [rad]
/// - slant_range = slant range [m]
///
/// ## Outputs:
/// - n = north coordinate [m] of input location from reference geodetic location
/// - e = east coordinate [m] of input location from reference geodetic location
/// - d = down coordinate [m] of input location from reference geodetic location
pub fn aer2ned(az: f64, el: f64, slant_range: f64) -> (f64, f64, f64) {
    let enu = aer2enu(az, el, slant_range);
    (enu.1, enu.0, -enu.2)
}

/// Returns the tuple (latitude,longitude,altitude) of coordinates in the Geodetic system
///
/// ## Inputs:
/// - az = azimuth angle [rad]
/// - el = elevation angle [rad]
/// - slant_range = slant range [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - lat = input latitude [rad]
/// - lon = input longitude [rad]
/// - alt = input altitude [m]
pub fn aer2geodetic(
    az: f64,
    el: f64,
    slant_range: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x, y, z) = aer2ecef(az, el, slant_range, lat0, lon0, alt0, r_ellips);
    ecef2geodetic(x, y, z, r_ellips)
}

/// Returns the tuple (u,v,w) of coordinates in the local vector system
///
/// ## Inputs:
/// - e = east coordinate [m] from reference geodetic location
/// - n = north coordinate [m] from reference geodetic location
/// - u = up coordinate [m] from reference geodetic location
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
///
/// ## Outputs:
/// - u = tangent vector component
/// - v = tangent vector component
/// - w = tangent vector component
pub fn enu2uvw(e: f64, n: f64, up: f64, lat0: f64, lon0: f64) -> (f64, f64, f64) {
    let t = lat0.cos() * up - lat0.sin() * n;

    let u = lon0.cos() * t - lon0.sin() * e;
    let v = lon0.sin() * t + lon0.cos() * e;
    let w = lat0.sin() * up + lat0.cos() * n;
    (u, v, w)
}

/// Returns the tuple (azimuth,elevation,slant range) of coordinates in the AER system
///
/// ## Inputs:
/// - e = east coordinate [m] from reference geodetic location
/// - n = north coordinate [m] from reference geodetic location
/// - u = up coordinate [m] from reference geodetic location
///
/// ## Outputs:
/// - az = azimuth angle [rad] of input location from reference geodetic location
/// - el = elevation angle [rad] of input location from reference geodetic location
/// - slant_range = slant range [m] of input location from reference geodetic location
pub fn enu2aer(e: f64, n: f64, u: f64) -> (f64, f64, f64) {
    let r = (e * e + n * n).sqrt();

    let slant_range = (r * r + u * u).sqrt();
    let el = u.atan2(r);
    let az = e.atan2(n).rem_euclid(2.0 * std::f64::consts::PI);

    (az, el, slant_range)
}

/// Returns the tuple (x,y,z) of coordinates in the ECEF system
///
/// ## Inputs:
/// - e = east coordinate [m] from reference geodetic location
/// - n = north coordinate [m] from reference geodetic location
/// - u = up coordinate [m] from reference geodetic location
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
pub fn enu2ecef(
    e: f64,
    n: f64,
    u: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x0, y0, z0) = geodetic2ecef(lat0, lon0, alt0, r_ellips);
    let (dx, dy, dz) = enu2uvw(e, n, u, lat0, lon0);

    (x0 + dx, y0 + dy, z0 + dz)
}

/// Returns the tuple (latitude,longitude,altitude) of coordinates in the Geodetic system
///
/// ## Inputs:
/// - e = east coordinate [m]  from reference geodetic location
/// - n = north coordinate [m]  from reference geodetic location
/// - u = up coordinate [m]  from reference geodetic location
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
/// ## Outputs:
/// - lat = latitude [rad]
/// - lon = longitude [rad]
/// - alt = altitude [m]
pub fn enu2geodetic(
    e: f64,
    n: f64,
    u: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x, y, z) = enu2ecef(e, n, u, lat0, lon0, alt0, r_ellips);
    let (lat, lon, alt) = ecef2geodetic(x, y, z, r_ellips);

    (lat, lon, alt)
}

/// Returns the tuple (x,y,z) of coordinates in the ECI system
///
/// ## Inputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
///
/// ## Outputs:
/// - x = x ECI coordinate [m]
/// - y = y ECI coordinate [m]
/// - z = z ECI coordinate [m]
pub fn ecef2eci(gst: f64, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let arr = matmul3(transpose3(r3(gst)), [x, y, z]);
    (arr[0], arr[1], arr[2])
}

/// Returns the tuple (latitude,longitude,altitude) of coordinates in the Geodetic system
/// with desired reference frame. Returns NaN for lat, lon, alt if the input ECEF coordinates are at the center of the Earth.
///
/// ## Inputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
/// - r_ellips = reference ellipsoid, defaults to WGS84
///
/// ## Outputs:
/// - lat = latitude [rad]
/// - lon = longitude [rad]
/// - alt = altitude [m]
pub fn ecef2geodetic(x: f64, y: f64, z: f64, r_ellips: Ellipsoid) -> (f64, f64, f64) {
    let (major, minor, _, _) = r_ellips.parameters();

    let q = (x * x + y * y).sqrt();

    // Handle poles case where q is zero to avoid division by zero
    if q <= f64::EPSILON {
        if z > 0.0 {
            return (std::f64::consts::FRAC_PI_2, 0.0, z - minor);
        }
        if z < 0.0 {
            return (-std::f64::consts::FRAC_PI_2, 0.0, -z - minor);
        }
        return (f64::NAN, f64::NAN, f64::NAN); // Indeterminate case at the center of the Earth
    }

    let r = (x * x + y * y + z * z).sqrt();
    let e = (major * major - minor * minor).sqrt();
    let var = r * r - e * e;
    let u = (0.5 * var + 0.5 * (var * var + 4.0 * e * e * z * z).sqrt()).sqrt();

    let hu_e = (u * u + e * e).sqrt();

    // FIX: this should improve stability near poles.
    let mut beta = (hu_e * z).atan2(u * q);
    let cos_beta = beta.cos();
    if cos_beta.abs() > 1e-12 {
        let sin_beta = beta.sin();
        let eps = ((minor * u - major * hu_e + e * e) * sin_beta)
            / (major * hu_e / cos_beta - e * e * cos_beta);
        beta += eps;
    }

    let lat = (major / minor * beta.tan()).atan();
    let lon = y.atan2(x);

    let v1 = z - minor * beta.sin();
    let v2 = q - major * beta.cos();

    let inside = (x * x / major / major) + (y * y / major / major) + (z * z / minor / minor) < 1.0;
    let alt = if inside {
        -(v1 * v1 + v2 * v2).sqrt()
    } else {
        (v1 * v1 + v2 * v2).sqrt()
    };

    (lat, lon, alt)
}

/// Returns the tuple (east,north,up) of coordinates in the ENU system
///
/// ## Inputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - e = east coordinate [m] of input ECEF location from reference geodetic location
/// - n = north coordinate [m] of input ECEF location from reference geodetic location
/// - u = up coordinate [m] of input ECEF location from reference geodetic location
pub fn ecef2enu(
    x: f64,
    y: f64,
    z: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x0, y0, z0) = geodetic2ecef(lat0, lon0, alt0, r_ellips);
    let (e, n, u) = uvw2enu(x - x0, y - y0, z - z0, lat0, lon0);
    (e, n, u)
}

/// Returns the tuple (north,east,down) of coordinates in the NED system
///
/// ## Inputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - n = north coordinate [m] of input location from reference geodetic location
/// - e = east coordinate [m] of input location from reference geodetic location
/// - d = down coordinate [m] of input location from reference geodetic location
pub fn ecef2ned(
    x: f64,
    y: f64,
    z: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let enu = ecef2enu(x, y, z, lat0, lon0, alt0, r_ellips);
    (enu.1, enu.0, -enu.2)
}

/// Returns the tuple (east,north,up) of coordinates in the ENU system
///
/// ## Inputs:
/// - u = tangent vector component
/// - v = tangent vector component
/// - w = tangent vector component
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
///
/// ## Outputs:
/// - e = east coordinate [m] of input location from reference geodetic location
/// - n = north coordinate [m] of input location from reference geodetic location
/// - u = up coordinate [m] of input location from reference geodetic location
pub fn uvw2enu(u: f64, v: f64, w: f64, lat0: f64, lon0: f64) -> (f64, f64, f64) {
    let t = lon0.cos() * u + lon0.sin() * v;
    let e = -lon0.sin() * u + lon0.cos() * v;
    let n = -lat0.sin() * t + lat0.cos() * w;
    let up = lat0.cos() * t + lat0.sin() * w;
    (e, n, up)
}

/// Returns the tuple (azimuth,elevation,slant range) of coordinates in the AER system
///
/// ## Inputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - az = azimuth angle [rad] of input location from reference geodetic location
/// - el = elevation angle [rad] of input location from reference geodetic location
/// - slant_range = slant range [m] of input location from reference geodetic location
pub fn ecef2aer(
    x: f64,
    y: f64,
    z: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (e, n, u) = ecef2enu(x, y, z, lat0, lon0, alt0, r_ellips);
    let (az, el, slant_range) = enu2aer(e, n, u);

    (az, el, slant_range)
}

/// Returns the tuple (azimuth,elevation,slant range) of coordinates in the AER system
///  
/// ## Inputs:
/// - x = x ECI coordinate [m]
/// - y = y ECI coordinate [m]
/// - z = z ECI coordinate [m]
/// - lat = reference latitude [rad]
/// - lon = reference longitude [rad]
/// - alt = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - az = azimuth angle [rad] of input location from reference geodetic location
/// - el = elevation angle [rad] of input location from reference geodetic location
/// - slant_range = slant range [m] of input location from reference geodetic location
pub fn eci2aer(
    gst: f64,
    x: f64,
    y: f64,
    z: f64,
    lat: f64,
    lon: f64,
    alt: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    let (x, y, z) = eci2ecef(gst, x, y, z);
    let (az, el, slant_range) = ecef2aer(x, y, z, lat, lon, alt, r_ellips);
    (az, el, slant_range)
}

/// Returns the tuple (x,y,z) of coordinates in the ECEF system
///
/// ## Inputs:
/// - gst = greenwich sidereal time
/// - x = x ECI coordinate [m]
/// - y = y ECI coordinate [m]
/// - z = z ECI coordinate [m]
///
/// ## Outputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
pub fn eci2ecef(gst: f64, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let arr = matmul3(r3(gst), [x, y, z]);
    (arr[0], arr[1], arr[2])
}

/// Returns the tuple (azimuth,elevation,slant range) of coordinates in the AER system
///
/// ## Inputs:
/// - n = north coordinate [m] of input location from reference geodetic location
/// - e = east coordinate [m] of input location from reference geodetic location
/// - d = down coordinate [m] of input location from reference geodetic location
///
/// ## Outputs:
/// - az = azimuth angle [rad] of input location from reference geodetic location
/// - el = elevation angle [rad] of input location from reference geodetic location
/// - slant_range = slant range [m] of input location from reference geodetic location
pub fn ned2aer(n: f64, e: f64, d: f64) -> (f64, f64, f64) {
    enu2aer(e, n, -d)
}

/// Returns the tuple (latitude,longitude,altitude) of coordinates in the Geodetic system
///
/// ## Inputs:
/// - n = north coordinate [m] of input location from reference geodetic location
/// - e = east coordinate [m] of input location from reference geodetic location
/// - d = down coordinate [m] of input location from reference geodetic location
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - lat = latitude [rad]
/// - lon = longitude [rad]
/// - alt = altitude [m]
pub fn ned2geodetic(
    n: f64,
    e: f64,
    d: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    enu2geodetic(e, n, -d, lat0, lon0, alt0, r_ellips)
}

/// Returns the tuple (x,y,z) of coordinates in the ECEF system
///
/// ## Inputs:
/// - n = north coordinate [m] from reference geodetic location
/// - e = east coordinate [m] from reference geodetic location
/// - d = down coordinate [m] from reference geodetic location
/// - lat0 = reference latitude [rad]
/// - lon0 = reference longitude [rad]
/// - alt0 = reference altitude [m]
/// - r_ellips = reference ellipsoid, defaults ref. frame is WGS84
///
/// ## Outputs:
/// - x = x ECEF coordinate [m]
/// - y = y ECEF coordinate [m]
/// - z = z ECEF coordinate [m]
pub fn ned2ecef(
    n: f64,
    e: f64,
    d: f64,
    lat0: f64,
    lon0: f64,
    alt0: f64,
    r_ellips: Ellipsoid,
) -> (f64, f64, f64) {
    enu2ecef(e, n, -d, lat0, lon0, alt0, r_ellips)
}

/// Returns the array result of 3-by-3-matrix that multiplies a 3-by-1 column array
pub fn matmul3(matrix: [f64; 9], col: [f64; 3]) -> [f64; 3] {
    [
        matrix[0] * col[0] + matrix[1] * col[1] + matrix[2] * col[2],
        matrix[3] * col[0] + matrix[4] * col[1] + matrix[5] * col[2],
        matrix[6] * col[0] + matrix[7] * col[1] + matrix[8] * col[2],
    ]
}

/// Returns the array representing a 3-by-3 rotation matrix of the input
pub fn r3(x: f64) -> [f64; 9] {
    [x.cos(), x.sin(), 0.0, -x.sin(), x.cos(), 0.0, 0.0, 0.0, 1.0]
}

/// Returns the array representing the transpose of the input 3-by-3 matrix
pub fn transpose3(x: [f64; 9]) -> [f64; 9] {
    [x[0], x[3], x[6], x[1], x[4], x[7], x[2], x[5], x[8]]
}

#[derive(Debug, Copy, Clone)]
pub enum Ellipsoid {
    /// WGS84: GPS Ellipsoid frame  
    /// semi-major axis: 6378137.0 [m]  
    /// flattening: 1.0/298.2572235630
    WGS84,
    /// WGS72: semi-major axis: 6378135.0 [m]    
    /// flattening: 1.0/298.26
    WGS72,
    /// WGS66: semi-major axis: 6378145.0 [m]    
    /// flattening: 1.0/298.25
    WGS66,
    /// WGS60: semi-major axis: 6378165.0 [m]    
    /// flattening: 1.0/298.3
    WGS60,
    /// PZ90: Glonass Ellipsoid frame  
    /// semi-major axis: 6378136.0 [m]
    /// flattening: 1/298.257839303
    PZ90,
    /// BDC, also known as CGCS2000,
    /// is the reference frame used by the
    /// Beidou constellation.  
    /// Semi-major axis: 6378137.0 [m]
    /// flattening: 1/298.257222101
    BDC,
    /// GRS80 reference ellipsoid  
    /// semi-major axis: 6378137.0 [m]  
    /// flattening: 1.0/298.257222101
    GRS80,
    /// Bessel reference ellipsoid   
    /// semi-major axis: 6377397.155 [m]
    /// flattening: 1.0/299.1528128
    Bessel,
    /// Airy reference ellipsoid   
    /// semi-major axis: 6377563.396 [m]  
    /// flattening: 1.0/299.3249646
    Airy,
    /// International reference ellipsoid   
    /// semi-major axis: 6378388.0 [m]  
    /// flattening: 1.0/297.0
    International,
}

impl Default for Ellipsoid {
    fn default() -> Ellipsoid {
        Ellipsoid::WGS84
    }
}

impl Ellipsoid {
    /// Returns the tuple representing the Ellipsoid frame.
    ///
    /// ## Outputs:
    /// - tuple.0 = semi-major axis [m]
    /// - tuple.1 = semi-minor axis [m]
    /// - tuple.2 = flattening [-]
    /// - tuple.3 = squared eccentricity [rad^2]
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        let (major, flattening): (f64, f64) = match self {
            Ellipsoid::WGS84 => (6378137.0, 1.0 / 298.257223563),
            Ellipsoid::WGS72 => (6378135.0, 1.0 / 298.26),
            Ellipsoid::WGS66 => (6378145.0, 1.0 / 298.25),
            Ellipsoid::WGS60 => (6378165.0, 1.0 / 298.3),
            Ellipsoid::PZ90 => (6378136.0, 1.0 / 298.257839303),
            Ellipsoid::BDC => (6378137.0, 1.0 / 298.257222101),
            Ellipsoid::GRS80 => (6378137.0, 1.0 / 298.2572221009),
            Ellipsoid::Bessel => (6377397.155, 1.0 / 299.1528128),
            Ellipsoid::Airy => (6377563.396, 1.0 / 299.3249646),
            Ellipsoid::International => (6378388.0, 1.0 / 297.0),
        };

        let minor = major * (1.0 - flattening);
        let ecc_sq = ((major * major) - (minor * minor)) / (major * major);
        (major, minor, flattening, ecc_sq)
    }
}

/// Returns the normal radius based on given latitude
/// and desired reference frame
pub fn get_radius_normal(lat: f64, r_ellips: Ellipsoid) -> f64 {
    let (major, _, _, squared_eccentricity) = r_ellips.parameters();
    major / ((1.0 - squared_eccentricity * lat.sin() * lat.sin()).sqrt())
}

/// Returns the radians [rad] value of the decimal degree [deg] input
#[deprecated(
    since = "0.1.6",
    note = "conversion from degrees to radians is natively supported. Use value.to_radians(). See https://doc.rust-lang.org/std/primitive.f64.html#method.to_radians"
)]
pub fn deg2rad(x: f64) -> f64 {
    x.to_radians()
}

/// Returns the decimal degree [deg] value of the radians [rad] input
#[deprecated(
    since = "0.1.6",
    note = "conversion from radians to degrees is natively supported. Use value.to_degrees(). See https://doc.rust-lang.org/std/primitive.f64.html#method.to_degrees"
)]
pub fn rad2deg(x: f64) -> f64 {
    x.to_degrees()
}

/// Returns the GST time as f64
///
/// ## Input
/// UTC time defined as: [year,month,day,hour,minute,second]
///
/// ## Output
/// Gst time as f64
pub fn utc2gst(utc: [i32; 6]) -> f64 {
    let mut year = utc[0] as f64;
    let mut month = utc[1] as f64;
    let day = utc[2] as f64;
    let h = utc[3] as f64;
    let m = utc[4] as f64;
    let s = utc[5] as f64;

    if month < 3.0 {
        year -= 1.0;
        month += 12.0;
    }

    let a = (year / 100.0).trunc();

    let b = 2.0 - a + (a / 4.0).trunc();

    let c = ((s / 60.0 + m) / 60.0 + h) / 24.0;

    let jd = (365.25 * (year + 4716.0)).trunc() + (30.6001 * (month + 1.0)).trunc() + day + b
        - 1524.5
        + c;

    let t_ut1 = (jd - 2451545.0) / 36525.0;

    let gmst_sec = 67310.54841 + 3.164400184812866e+09 * t_ut1 + 0.093104 * t_ut1 * t_ut1
        - 6.2e-6 * t_ut1 * t_ut1 * t_ut1;

    (gmst_sec * 2.0 * std::f64::consts::PI / 86400.0).rem_euclid(2.0 * std::f64::consts::PI)
}

/// Return the round toward zero value of the input
#[deprecated(
    since = "0.1.7",
    note = "rounding toward zero is natively supported. Use value.trunc(). See https://doc.rust-lang.org/std/primitive.f64.html#method.trunc"
)]
pub fn fix(x: f64) -> f64 {
    let mut out = x;
    if out < 0.0 {
        out = x.ceil();
    } else {
        out = x.floor();
    }
    out
}

/// Earth radius (m)
pub const EARTH_RADIUS: f64 = 6371E3_f64;

/// Returns distance (m) between two decimal degrees coordinates::
/// coord1: (lat,lon), coord2: (lat, lon)
pub fn distance(coord1: (f64, f64), coord2: (f64, f64)) -> f64 {
    let dphi = coord2.0.to_radians() - coord1.0.to_radians();
    let d_lambda = coord2.1.to_radians() - coord1.1.to_radians();
    let a: f64 = (dphi / 2.0_f64).sin().powf(2.0_f64)
        + coord1.0.to_radians().cos()
            * coord2.0.to_radians().cos()
            * (d_lambda / 2.0_f64).sin().powf(2.0_f64);

    let a = a.clamp(0.0, 1.0); // Clamp a to the range [0, 1] to prevent numerical issues
    let c = 2.0_f64 * a.sqrt().atan2((1.0 - a).sqrt());
    EARTH_RADIUS * c
}

/// Returns the string geohash of the geodetic coordinates (latitude, longitude) with desired precision.
/// Precision is clamped to the range [0, 12].
/// ## Inputs:
/// - lat = latitude [rad]
/// - lon = longitude [rad]
/// - precision = number of characters in the geohash string, defaults to 12
///
/// ## Outputs:
/// - geohash string of the input geodetic coordinates with desired precision
pub fn geodetic2geohash(lat: f64, lon: f64, precision: usize) -> String {
    // Geohash typically uses degrees for latitude and longitude
    let lat = lat.to_degrees();
    let lon = lon.to_degrees();
    let precision = precision.clamp(0, 12);

    let mut lat_range = (-90.0, 90.0);
    let mut lon_range = (-180.0, 180.0);
    let mut bits = 0u8; // 5 bits per character
    let mut hash = String::with_capacity(precision);

    let charset = b"0123456789bcdefghjkmnpqrstuvwxyz";

    for i in 0..(precision * 5) {
        let is_even = i % 2 == 0;
        let (val, range) = if is_even {
            (lon, &mut lon_range)
        } else {
            (lat, &mut lat_range)
        };

        let mid = (range.0 + range.1) / 2.0;
        bits <<= 1;
        if val >= mid {
            bits |= 1;
            range.0 = mid;
        } else {
            range.1 = mid;
        }

        if (i + 1) % 5 == 0 {
            hash.push(charset[bits as usize] as char);
            bits = 0;
        }
    }
    hash
}

// Module for unit tests
#[cfg(test)]
mod tests;
