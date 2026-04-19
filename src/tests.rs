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
#[allow(deprecated)]
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

    let (x, y, z) = aer2eci(
        gst,
        az,
        el,
        slant_range,
        lat0,
        lon0,
        alt0,
        Ellipsoid::default(),
    );

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

    let (lat, lon, alt) = aer2geodetic(az, el, slant_range, lat0, lon0, alt0, Ellipsoid::default());

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

    let (x, y, z) = aer2ecef(
        azref,
        elref,
        rangeref,
        lat0,
        lon0,
        alt0,
        Ellipsoid::default(),
    );
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
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::WGS66.parameters();
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::WGS60.parameters();
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::PZ90.parameters();
    assert!((a - 6378136.0).abs() < 1E-6);
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((1.0 / f - 298.257839303).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::GRS80.parameters();
    assert!((a - 6378137.0).abs() < 1E-6);
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((1.0 / f - 298.257222101).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::BDC.parameters();
    assert!((a - 6378137.0).abs() < 1E-6);
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((1.0 / f - 298.257222101).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::Bessel.parameters();
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::International.parameters();
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
    let (a, b, f, e) = Ellipsoid::Airy.parameters();
    assert!((b - a * (1.0 - f)).abs() < 1E-6);
    assert!((e - (f * (2.0 - f))).abs() < 1E-6);
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

#[test]
fn test_ecef2geodetic_poles() {
    let (_, minor, _, _) = Ellipsoid::default().parameters();

    let (lat, lon, alt) = ecef2geodetic(0.0, 0.0, minor, Ellipsoid::default());
    assert!((lat - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    assert!(lon.abs() < 1e-12);
    assert!(alt.abs() < 1e-12);

    let (lat, lon, alt) = ecef2geodetic(0.0, 0.0, -minor, Ellipsoid::default());
    assert!((lat + std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    assert!(lon.abs() < 1e-12);
    assert!(alt.abs() < 1e-12);
}

#[test]
fn test_ecef2geodetic_near_pole() {
    let lat_ref = 89.999999_f64.to_radians();
    let lon_ref = 30.0_f64.to_radians();
    let alt_ref = 1000.0;

    let (x, y, z) = geodetic2ecef(lat_ref, lon_ref, alt_ref, Ellipsoid::default());
    let (lat, lon, alt) = ecef2geodetic(x, y, z, Ellipsoid::default());

    assert!((lat - lat_ref).abs() < 1e-8);
    assert!((lon - lon_ref).abs() < 1e-8);
    assert!((alt - alt_ref).abs() < 1e-5);
}

#[test]
fn test_distance_antipodal() {
    let d = distance((0.0, 0.0), (0.0, 180.0));
    assert!(d.is_finite());
}

#[test]
fn test_geohash_known_points() {
    let cities: Vec<(&str, f64, f64, usize, &str)> = vec![
        ("Tokyo", 35.6897, 139.6922, 5, "xn774"),
        ("London", 51.5074, -0.1278, 6, "gcpvj0"),
        ("Washington DC", 38.8951, -77.0364, 7, "dqcjqbx"),
        ("Brasília", -15.7975, -47.8919, 4, "6vjy"),
        ("Canberra", -35.2809, 149.1300, 8, "r3dp3931"),
        ("Seoul", 37.5665, 126.9780, 12, "wydm9qy89z5m"), //clamped to 12
    ];

    for (name, lat, lon, len, expected) in cities {
        let result = geodetic2geohash(lat.to_radians(), lon.to_radians(), len);
        assert_eq!(
            result, expected,
            "City {} failed. Expected {}, got {}",
            name, expected, result
        );
        assert_eq!(result.len(), len, "Length mismatch for {}", name);
    }
}

#[test]
fn test_geohash_boundary_conditions() {
    // Exact center (0,0) - Prime Meridian & Equator
    // Degrees: (0.0, 0.0) -> "s000..."
    let h_center = geodetic2geohash(0.0, 0.0, 1);
    assert_eq!(h_center, "s");

    // Exact North Pole
    let h_pole = geodetic2geohash(std::f64::consts::PI / 2.0, 0.0, 5);
    assert!(h_pole.len() == 5);
}
