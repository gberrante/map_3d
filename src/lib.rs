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
pub fn geodetic2aer(lat: f64, lon: f64, alt: f64, lat0: f64, lon0: f64, alt0: f64,
                    r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn geodetic2enu(lat: f64, lon: f64, alt: f64, lat0: f64, lon0: f64, alt0: f64,
                    r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn geodetic2ned(lat: f64, lon: f64, alt: f64, lat0: f64, lon0: f64, alt0: f64,
                    r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn aer2ecef(az: f64, el: f64, slant_range: f64, lat0: f64, lon0: f64, alt0: f64,
                r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn aer2eci(gst: f64, az: f64, el: f64, slant_range: f64, lat0: f64, lon0: f64, alt0: f64,
               r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn aer2geodetic(az: f64, el: f64, slant_range: f64, lat0: f64, lon0: f64, alt0: f64,
                    r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn enu2uvw(et: f64, nt: f64, up: f64, lat0: f64, lon0: f64) -> (f64, f64, f64) {
    let t = lat0.cos() * up - lat0.sin() * nt;

    let u = lon0.cos() * t - lon0.sin() * et;
    let v = lon0.sin() * t + lon0.cos() * et;
    let w = lat0.sin() * up + lat0.cos() * nt;
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
pub fn enu2ecef(e: f64, n: f64, u: f64, lat0: f64, lon0: f64, alt0: f64,
                r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn enu2geodetic(e: f64, n: f64, u: f64, lat0: f64, lon0: f64, alt0: f64,
                    r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
/// with desired reference frame.
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

    let r = (x * x + y * y + z * z).sqrt();
    let e = (major * major - minor * minor).sqrt();
    let var = r * r - e * e;
    let u = (0.5 * var + 0.5 * (var * var + 4.0 * e * e * z * z).sqrt()).sqrt();

    let q = (x * x + y * y).sqrt();
    let hu_e = (u * u + e * e).sqrt();
    let mut beta = (hu_e / u * z / q).atan();

    let eps = ((minor * u - major * hu_e + e * e) * beta.sin())
        / (major * hu_e / beta.cos() - e * e * beta.cos());
    beta += eps;

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
pub fn ecef2enu(x: f64, y: f64, z: f64, lat0: f64, lon0: f64, alt0: f64,
                r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn ecef2ned(x: f64, y: f64, z: f64, lat0: f64, lon0: f64, alt0: f64,
                r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
    let u = lat0.cos() * t + lat0.sin() * w;
    (e, n, u)
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
pub fn ecef2aer(x: f64, y: f64, z: f64, lat0: f64, lon0: f64, alt0: f64,
                r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn eci2aer(gst: f64, x: f64, y: f64, z: f64, lat: f64, lon: f64, alt: f64,
               r_ellips: Ellipsoid) -> (f64, f64, f64) {
    let (x, y, z) = eci2ecef(gst, x, y, z);
    let (az, el, slant_range) = ecef2aer(x, y, z, lat, lon, alt, r_ellips);
    (az, el, slant_range)
}

/// Returns the tuple (x,y,z) of coordinates in the ECEF system
///
/// ## Inputs:
/// - gst = greenwhich sidereal time
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
pub fn ned2geodetic(n: f64, e: f64, d: f64, lat0: f64, lon0: f64, alt0: f64,
                    r_ellips: Ellipsoid) -> (f64, f64, f64) {
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
pub fn ned2ecef(n: f64, e: f64, d: f64, lat0: f64, lon0: f64, alt0: f64,
                r_ellips: Ellipsoid) -> (f64, f64, f64) {
    enu2ecef(e, n, -d, lat0, lon0, alt0, r_ellips)
}

/// Returns the array result of 3-by-3-matrix that multiplies a 3-by-1 column array
pub fn matmul3(matrix: [f64; 9], col: [f64; 3]) -> [f64; 3] {
    let out: [f64; 3] = [
        matrix[0] * col[0] + matrix[1] * col[1] + matrix[2] * col[2],
        matrix[3] * col[0] + matrix[4] * col[1] + matrix[5] * col[2],
        matrix[6] * col[0] + matrix[7] * col[1] + matrix[8] * col[2],
    ];
    out
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
            Ellipsoid::Bessel => (6377397.155, 299.1528128),
            Ellipsoid::Airy => (6377563.396, 299.3249646),
            Ellipsoid::International => (6378388.0, 297.0),
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
#[deprecated(since = "0.1.6", note = "conversion from degrees to radians is natively supported. Use value.to_radians(). See https://doc.rust-lang.org/std/primitive.f64.html#method.to_radians")]
pub fn deg2rad(x: f64) -> f64 {
    x.to_radians()
}

/// Returns the decimal degree [deg] value of the radians [rad] input
#[deprecated(since = "0.1.6", note = "conversion from radians to degrees is natively supported. Use value.to_degrees(). See https://doc.rust-lang.org/std/primitive.f64.html#method.to_degrees")]
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

    let a = fix(year / 100.0);

    let b = 2.0 - a + fix(a / 4.0);

    let c = ((s / 60.0 + m) / 60.0 + h) / 24.0;

    let jd = fix(365.25 * (year + 4716.0)) + fix(30.6001 * (month + 1.0)) + day + b - 1524.5 + c;

    let t_ut1 = (jd - 2451545.0) / 36525.0;

    let gmst_sec = 67310.54841 + 3.164400184812866e+09 * t_ut1 + 0.093104 * t_ut1 * t_ut1
        - 6.2e-6 * t_ut1 * t_ut1 * t_ut1;

    (gmst_sec * 2.0 * std::f64::consts::PI / 86400.0).rem_euclid(2.0 * std::f64::consts::PI)
}

/// Return the round toward zero value of the input
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
    let c = 2.0_f64 * a.powf(0.5_f64).atan2((1.0 - a).powf(0.5_f64));
    EARTH_RADIUS * c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_utc2gst() {
        let datetime: [i32; 6] = [2020, 5, 12, 18, 2, 10];
        let gst_ref = 2.469809475597415;

        let gst = utc2gst(datetime);

        assert!((gst - gst_ref).abs() < 1e-8);

        let datetime2: [i32; 6] = [2020, 1, 12, 18, 2, 10];
        let gst_ref2 = 0.388271658105431;

        let gst2 = utc2gst(datetime2);

        assert!((gst2 - gst_ref2).abs() < 1e-8);
    }

    #[test]
    fn test_fix() {
        let x1 = 3.7;
        let x2 = -4.67;

        assert_eq!(fix(x1), 3.0);
        assert_eq!(fix(x2), -4.0);
    }

    #[test]
    fn test_geodetic2ecef() {
        let lat = 30.14988205_f64.to_radians();
        let lon = 91.38733072_f64.to_radians();
        let alt = 4031.0;

        let (x, y, z) = geodetic2ecef(lat, lon, alt, Ellipsoid::default());

        let xref = -1.337281037300386e+05;
        let yref = 5.521796910920261e+06;
        let zref = 3.186776473672415e+06;

        assert!((x - xref).abs() < 1e-3);
        assert!((y - yref).abs() < 1e-3);
        assert!((z - zref).abs() < 1e-3);
    }

    #[test]
    fn test_geodetic2aer() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;

        let lat = 42.002581974253744_f64.to_radians();
        let lon = -81.997751960067460_f64.to_radians();
        let alt = 1.139701799575106e+03;

        let azref = 32.999999999989740_f64.to_radians();
        let elref = 69.999999999945540_f64.to_radians();
        let rangeref = 1000.0;

        let (a, e, r) = geodetic2aer(lat, lon, alt, lat0, lon0, alt0, Ellipsoid::default());

        assert!((a - azref).abs() < 1e-3);
        assert!((e - elref).abs() < 1e-3);
        assert!((r - rangeref).abs() < 1e-3);
    }

    #[test]
    fn test_geodetic2enu() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;

        let lat = 42.002581974253744_f64.to_radians();
        let lon = -81.997751960067460_f64.to_radians();
        let alt = 1.139701799575106e+03;

        let eref = 1.862775208168244e+02;
        let nref = 2.868422278521820e+02;
        let uref = 9.396926207845534e+02;

        let (e, n, u) = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, Ellipsoid::default());

        assert!((e - eref).abs() < 1e-3);
        assert!((n - nref).abs() < 1e-3);
        assert!((u - uref).abs() < 1e-3);
    }

    #[test]
    fn test_aer2ecef() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;

        let az = 33.0_f64.to_radians();
        let el = 70.0_f64.to_radians();
        let slant_range = 1000.0;

        let (x, y, z) = aer2ecef(az, el, slant_range, lat0, lon0, alt0, Ellipsoid::default());

        let xref = 6.609301927610816e+05;
        let yref = -4.701424222957011e+06;
        let zref = 4.246579604632881e+06;

        assert!((x - xref).abs() < 1e-3);
        assert!((y - yref).abs() < 1e-3);
        assert!((z - zref).abs() < 1e-3);
    }

    #[test]
    fn test_aer2enu() {
        let az = 33.0_f64.to_radians();
        let el = 70.0_f64.to_radians();
        let slant_range = 1000.0;

        let eref = 1.862775208165935e+02;
        let nref = 2.868422278517140e+02;
        let uref = 9.396926207859083e+02;

        let (e, n, u) = aer2enu(az, el, slant_range);

        assert!((e - eref).abs() < 1e-3);
        assert!((n - nref).abs() < 1e-3);
        assert!((u - uref).abs() < 1e-3);
    }

    #[test]
    fn test_aer2eci() {
        let az = 162.55_f64.to_radians();
        let el = 55.12_f64.to_radians();
        let slant_range = 384013940.9;
        let gst = 4.501012562811752;

        let lat0 = 28.4_f64.to_radians();
        let lon0 = -80.5_f64.to_radians();
        let alt0 = 2.7;

        let xref = -3.849714979138141e+08;
        let yref = -4.836588977863766e+07;
        let zref = -3.143285462295778e+07;

        let (x, y, z) = aer2eci(gst, az, el, slant_range, lat0, lon0, alt0, Ellipsoid::default());

        assert!((x - xref).abs() < 1e-3);
        assert!((y - yref).abs() < 1e-3);
        assert!((z - zref).abs() < 1e-3);
    }

    #[test]
    fn test_aer2geodetic() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;

        let az = 32.999999999989740_f64.to_radians();
        let el = 69.999999999945540_f64.to_radians();
        let slant_range = 1000.0;

        let latref = 42.002581974253744_f64.to_radians();
        let lonref = -81.997751960067460_f64.to_radians();
        let altref = 1.139701799575106e+03;

        let (lat, lon, alt) =
            aer2geodetic(az, el, slant_range, lat0, lon0, alt0, Ellipsoid::default());

        assert!((lat - latref).abs() < 1e-8);
        assert!((lon - lonref).abs() < 1e-8);
        assert!((alt - altref).abs() < 1e-8);
    }

    #[test]
    fn test_enu2aer() {
        let e = 1.862775210000000e+02;
        let n = 2.868422200000000e+02;
        let u = 9.396926200000000e+02;

        let azref = 33.0_f64.to_radians();
        let elref = 70.0_f64.to_radians();
        let rangeref = 1000.0;

        let (az, el, range) = enu2aer(e, n, u);

        assert!((az - azref).abs() < 1e-3);
        assert!((el - elref).abs() < 1e-3);
        assert!((range - rangeref).abs() < 1e-3);
    }

    #[test]
    fn test_enu2ecef() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;
        let e = 1.862775210000000e+02;
        let n = 2.868422200000000e+02;
        let u = 9.396926200000000e+02;

        let xref = 6.609301927610815e+05;
        let yref = -4.701424222957011e+06;
        let zref = 4.246579604632881e+06;

        let (x, y, z) = enu2ecef(e, n, u, lat0, lon0, alt0, Ellipsoid::default());

        assert!((x - xref).abs() < 1e-3);
        assert!((y - yref).abs() < 1e-3);
        assert!((z - zref).abs() < 1e-3);
    }

    #[test]
    fn test_enu2geodetic() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;
        let e = 0.0;
        let n = 0.0;
        let u = -1.0;

        let latref = 41.999999999999993_f64.to_radians();
        let lonref = -82.0_f64.to_radians();
        let altref = 1.990000000007368e+02;

        let (lat, lon, alt) = enu2geodetic(e, n, u, lat0, lon0, alt0, Ellipsoid::default());

        assert!((lat - latref).abs() < 1e-8);
        assert!((lon - lonref).abs() < 1e-8);
        assert!((alt - altref).abs() < 1e-8);
    }

    #[test]
    fn test_ecef2geodetic() {
        let latref = 30.14988205_f64.to_radians();
        let lonref = 91.38733072_f64.to_radians();
        let altref = 4031.0;

        let (x, y, z) = geodetic2ecef(latref, lonref, altref, Ellipsoid::default());
        let (lat, lon, alt) = ecef2geodetic(x, y, z, Ellipsoid::default());

        assert!((lat - latref).abs() < 1e-8);
        assert!((lon - lonref).abs() < 1e-8);
        assert!((alt - altref).abs() < 1e-8);

        let (x, y, z) = geodetic2ecef(latref, lonref, altref - 5000.0, Ellipsoid::default());
        let (lat, lon, alt) = ecef2geodetic(x, y, z, Ellipsoid::default());

        assert!((lat - latref).abs() < 1e-8);
        assert!((lon - lonref).abs() < 1e-8);
        assert!((alt - (altref - 5000.0)).abs() < 1e-8);
    }

    #[test]
    fn test_ecef2enu() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;
        let eref = 1.862775210000000e+02;
        let nref = 2.868422200000000e+02;
        let uref = 9.396926200000000e+02;

        let (x, y, z) = enu2ecef(eref, nref, uref, lat0, lon0, alt0, Ellipsoid::default());
        let (e, n, u) = ecef2enu(x, y, z, lat0, lon0, alt0, Ellipsoid::default());

        assert!((e - eref).abs() < 1e-3);
        assert!((n - nref).abs() < 1e-3);
        assert!((u - uref).abs() < 1e-3);
    }

    #[test]
    fn test_ecef2aer() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;

        let azref = 33.0_f64.to_radians();
        let elref = 70.0_f64.to_radians();
        let rangeref = 1000.0;

        let (x, y, z) = aer2ecef(azref, elref, rangeref, lat0, lon0, alt0, Ellipsoid::default());
        let (az, el, range) = ecef2aer(x, y, z, lat0, lon0, alt0, Ellipsoid::default());

        assert!((az - azref).abs() < 1e-3);
        assert!((el - elref).abs() < 1e-3);
        assert!((range - rangeref).abs() < 1e-3);
    }

    #[test]
    fn test_eci2aer() {
        let azref = 162.55_f64.to_radians();
        let elref = 55.12_f64.to_radians();
        let rangeref = 384013940.9;
        let gst = 4.501012562811752;

        let lat0 = 28.4_f64.to_radians();
        let lon0 = -80.5_f64.to_radians();
        let alt0 = 2.7;

        let (x, y, z) = aer2eci(
            gst,
            azref,
            elref,
            rangeref,
            lat0,
            lon0,
            alt0,
            Ellipsoid::default(),
        );
        let (az, el, range) = eci2aer(gst, x, y, z, lat0, lon0, alt0, Ellipsoid::default());

        assert!((az - azref).abs() < 1e-3);
        assert!((el - elref).abs() < 1e-3);
        assert!((range - rangeref).abs() < 1e-3);
    }

    #[test]
    fn test_ned2geodetic() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;
        let e = 0.0;
        let n = 0.0;
        let d = 1.0;

        let latref = 41.999999999999993_f64.to_radians();
        let lonref = -82.0_f64.to_radians();
        let altref = 1.990000000007368e+02;

        let (lat, lon, alt) = ned2geodetic(n, e, d, lat0, lon0, alt0, Ellipsoid::default());

        assert!((lat - latref).abs() < 1e-8);
        assert!((lon - lonref).abs() < 1e-8);
        assert!((alt - altref).abs() < 1e-8);
    }

    #[test]
    fn test_geodetic2ned() {
        let lat = 41.999999999999993_f64.to_radians();
        let lon = -82.0_f64.to_radians();
        let alt = 1.990000000007368e+02;
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;

        let eref = 0.0;
        let nref = 0.0;
        let dref = 1.0;

        let (n, e, d) = geodetic2ned(lat, lon, alt, lat0, lon0, alt0, Ellipsoid::default());

        assert!((e - eref).abs() < 1e-3);
        assert!((n - nref).abs() < 1e-3);
        assert!((d - dref).abs() < 1e-3);
    }

    #[test]
    fn test_aer2ned() {
        let az = 33.0_f64.to_radians();
        let el = 70.0_f64.to_radians();
        let slant_range = 1000.0;

        let eref = 1.862775208165935e+02;
        let nref = 2.868422278517140e+02;
        let dref = -9.396926207859083e+02;

        let (n, e, d) = aer2ned(az, el, slant_range);

        assert!((e - eref).abs() < 1e-3);
        assert!((n - nref).abs() < 1e-3);
        assert!((d - dref).abs() < 1e-3);
    }

    #[test]
    fn test_ned2aer() {
        let az_ref = 33.0_f64.to_radians();
        let el_ref = 70.0_f64.to_radians();
        let range_ref = 1000.0;

        let e = 1.862775208165935e+02;
        let n = 2.868422278517140e+02;
        let d = -9.396926207859083e+02;

        let (az, el, range) = ned2aer(n, e, d);

        assert!((az - az_ref).abs() < 1e-6);
        assert!((el - el_ref).abs() < 1e-6);
        assert!((range - range_ref).abs() < 1e-3);
    }

    #[test]
    fn test_ned2ecef() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;
        let e = 1.862775210000000e+02;
        let n = 2.868422200000000e+02;
        let d = -9.396926200000000e+02;

        let xref = 6.609301927610815e+05;
        let yref = -4.701424222957011e+06;
        let zref = 4.246579604632881e+06;

        let (x, y, z) = ned2ecef(n, e, d, lat0, lon0, alt0, Ellipsoid::default());

        assert!((x - xref).abs() < 1e-3);
        assert!((y - yref).abs() < 1e-3);
        assert!((z - zref).abs() < 1e-3);
    }

    #[test]
    fn test_ellipsoid_references() {
        let (a, b, f, e) = Ellipsoid::WGS84.parameters();
        assert!((a - 6378137.0).abs() < 1E-6);
        assert!((b - 6356752.314245).abs() < 1E-6);
        assert!((1.0 / f - 298.257223563).abs() < 1E-6);
        assert!((e - 6.6943799E-3_f64).abs() < 1E-6);
        let (a, b, f, e) = Ellipsoid::WGS72.parameters();
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::WGS66.parameters();
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::WGS60.parameters();
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::PZ90.parameters();
        assert!((a - 6378136.0).abs() < 1E-6);
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((1.0 / f - 298.257839303).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::GRS80.parameters();
        assert!((a - 6378137.0).abs() < 1E-6);
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((1.0 / f - 298.257222101).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::BDC.parameters();
        assert!((a - 6378137.0).abs() < 1E-6);
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((1.0 / f - 298.257222101).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::Bessel.parameters();
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::International.parameters();
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
        let (a, b, f, e) = Ellipsoid::Airy.parameters();
        assert!((b - a * (1.0 - f)).abs() < 1E-6);
        assert!((e - (f * (2.0 - f))) < 1E-6);
    }

    #[test]
    fn test_ecef2ned() {
        let lat0 = 42.0_f64.to_radians();
        let lon0 = -82.0_f64.to_radians();
        let alt0 = 200.0;
        let eref = 1.862775210000000e+02;
        let nref = 2.868422200000000e+02;
        let dref = -9.396926200000000e+02;

        let (x, y, z) = ned2ecef(nref, eref, dref, lat0, lon0, alt0, Ellipsoid::default());
        let (n, e, d) = ecef2ned(x, y, z, lat0, lon0, alt0, Ellipsoid::default());

        assert!((e - eref).abs() < 1e-3);
        assert!((n - nref).abs() < 1e-3);
        assert!((d - dref).abs() < 1e-3);
    }

    #[test]
    fn test_distance_calculation() {
        let new_york = (40.730610, -73.935242);
        let paris = (48.856614, 2.3522219);
        let buenos_aires = (-34.603722, -58.381592);
        let sydney = (-33.867487, 151.20699);
        // TEST 1
        let expected_km = 5831.0_f64;
        let d_km = distance(new_york, paris) / 1000.0_f64;
        assert!((expected_km - d_km).abs() < 1.0);
        // TEST2
        let expected_km = 8527.0_f64;
        let d_km = distance(new_york, buenos_aires) / 1000.0_f64;
        assert!((expected_km - d_km).abs() < 1.0);
        // TEST3
        let expected_km = 15990.0_f64;
        let d_km = distance(new_york, sydney) / 1000.0_f64;
        assert!((expected_km - d_km).abs() < 10.0);
        // TEST4
        let expected_km = 11050.0_f64;
        let d_km = distance(buenos_aires, paris) / 1000.0_f64;
        assert!((expected_km - d_km).abs() < 10.0);
    }
}
