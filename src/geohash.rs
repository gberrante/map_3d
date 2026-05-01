//! map_3d: geographic coordinates conversions

// Charset for geohash encoding (base32)
const GEOHASH_CHARSET: &[u8] = b"0123456789bcdefghjkmnpqrstuvwxyz";

/// Returns the string geohash of the geodetic coordinates (latitude, longitude) with desired precision.
/// Precision is clamped to the range [0, 12].
/// ## Inputs:
/// - lat = latitude [rad]
/// - lon = longitude [rad]
/// - precision = number of characters in the geohash string, defaults to 12
///
/// ## Outputs:
/// - geohash string of the input geodetic coordinates with desired precision
pub fn encode(lat: f64, lon: f64, precision: usize) -> String {
    // Geohash typically uses degrees for latitude and longitude
    let lat = lat.to_degrees();
    let lon = lon.to_degrees();
    let precision = precision.clamp(0, 12);

    let mut lat_range = (-90.0, 90.0);
    let mut lon_range = (-180.0, 180.0);
    let mut bits = 0u8; // 5 bits per character
    let mut hash = String::with_capacity(precision);

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
            hash.push(GEOHASH_CHARSET[bits as usize] as char);
            bits = 0;
        }
    }
    hash
}

/// Returns the geodetic coordinates (latitude, longitude) corresponding to the input geohash string.
/// ## Inputs:
/// - geohash = geohash string to decode
/// ## Outputs:
/// - (latitude, longitude, precision) corresponding to the input geohash string
/// Note: the returned latitude and longitude are the center of the geohash cell, not the exact coordinates.
pub fn decode(x: &str) -> Option<(f64, f64, usize)> {
    let mut lat_range = (-90.0, 90.0);
    let mut lon_range = (-180.0, 180.0);
    let mut is_even = true;

    for c in x.bytes() {
        let idx = GEOHASH_CHARSET.iter().position(|&x| x == c)?;
        for mask in [16, 8, 4, 2, 1].iter() {
            if is_even {
                if (idx & mask) != 0 {
                    lon_range.0 = (lon_range.0 + lon_range.1) / 2.0;
                } else {
                    lon_range.1 = (lon_range.0 + lon_range.1) / 2.0;
                }
            } else {
                if (idx & mask) != 0 {
                    lat_range.0 = (lat_range.0 + lat_range.1) / 2.0;
                } else {
                    lat_range.1 = (lat_range.0 + lat_range.1) / 2.0;
                }
            }
            is_even = !is_even;
        }
    }

    let lat: f64 = (lat_range.0 + lat_range.1) / 2.0;
    let lon: f64 = (lon_range.0 + lon_range.1) / 2.0;
    Some((lat.to_radians(), lon.to_radians(), x.len()))
}
