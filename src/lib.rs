pub fn geodetic2ecef(lat: f64, lon: f64, alt: f64) -> (f64,f64,f64){
    let n = get_radius_normal(lat);
    let (major,minor,_,_) = wgs84();

    let x = (n + alt) * lat.cos() * lon.cos();
    let y = (n + alt) * lat.cos() * lon.sin();
    let z = (n * (minor / major) * (minor / major) + alt) * lat.sin();

    (x,y,z)
}

pub fn geodetic2aer(lat: f64, lon: f64, alt: f64,lat0: f64, lon0: f64, alt0: f64) -> (f64,f64,f64){
    let (e,n,u) = geodetic2enu(lat,lon,alt,lat0,lon0,alt0);
    let (az,el,slant_range) = enu2aer(e, n, u);

    (az,el,slant_range)
}

pub fn geodetic2enu(lat: f64, lon: f64, alt: f64,lat0: f64, lon0: f64, alt0: f64) -> (f64,f64,f64){
    let (x1,y1,z1) = geodetic2ecef(lat, lon, alt);
    let (x2,y2,z2) = geodetic2ecef(lat0, lon0, alt0);

    let (e,n,u) = ecef2enuv(x1-x2, y1-y2, z1-z2, lat0, lon0);

    (e,n,u)
} 

// --------------------------------------

pub fn aer2ecef(az : f64, el: f64,slant_range :f64, lat0: f64, lon0: f64, alt0: f64)-> (f64,f64,f64){

    let (x0,y0,z0) = geodetic2ecef(lat0, lon0, alt0);
    let (e,n,u) = aer2enu(az,el,slant_range);
    let (dx,dy,dz) = enu2uvw(e,n,u,lat0,lon0);
    (x0+dx,y0+dy,z0+dz)
}

pub fn aer2enu(az : f64, el: f64,slant_range :f64) -> (f64,f64,f64){
    let r = slant_range*el.cos();
    (r*az.sin(),r*az.cos(),slant_range*el.sin())
}

pub fn aer2eci(gst :f64, az : f64, el: f64,slant_range :f64, lat0: f64, lon0: f64, alt0: f64) -> (f64,f64,f64){
    let (x1,y1,z1) = aer2ecef(az, el, slant_range, lat0, lon0, alt0);
    ecef2eci(gst, x1, y1, z1)
}

pub fn aer2geodetic(az : f64, el: f64,slant_range :f64, lat0: f64, lon0: f64, alt0: f64)-> (f64,f64,f64) {
    let (x,y,z) = aer2ecef(az, el, slant_range, lat0, lon0, alt0);
    ecef2geodetic(x,y,z)
}

// --------------------------------------

pub fn enu2uvw(et : f64, nt: f64,up :f64, lat0: f64, lon0: f64) -> (f64,f64,f64){
    let t = lat0.cos()*up - lat0.sin()*nt;

    let u = lon0.cos()*t-lon0.sin()*et;
    let v = lon0.sin()*t+lon0.cos()*et;
    let w = lat0.sin()*up+lat0.cos()*nt;
    (u,v,w)
}

pub fn enu2aer(e : f64, n: f64,u :f64) -> (f64,f64,f64){
    let r = (e*e+n*n).sqrt();

    let slant_range = (r*r+u*u).sqrt();
    let el = u.atan2(r);
    let az = e.atan2(n) % (2.0*std::f64::consts::PI);

    (az,el,slant_range)

}

pub fn enu2ecef(e : f64, n: f64,u :f64, lat0: f64, lon0: f64, alt0 : f64) -> (f64,f64,f64){
    let (x0,y0,z0) = geodetic2ecef(lat0, lon0, alt0);
    let (dx,dy,dz) = enu2uvw(e, n, u, lat0, lon0);

    (x0+dx,y0+dy,z0+dz)
}

pub fn enu2ecefv(e : f64, n: f64,u :f64, lat0: f64, lon0: f64) -> (f64,f64,f64){
    let t = lat0.cos() * u - lat0.sin() * n;
    let u = lon0.cos() * t - lon0.sin() * e;
    let v = lon0.sin() * t + lon0.cos() * e;
    let w = lat0.sin() * u + lat0.cos() * n;

    (u,v,w)
}

pub fn enu2geodetic(e : f64, n: f64,u :f64, lat0: f64, lon0: f64, alt0 : f64) -> (f64,f64,f64){
    let (x,y,z) = enu2ecef(e, n, u, lat0, lon0, alt0);
    let (lat,lon,alt) = ecef2geodetic(x, y, z);

    (lat,lon,alt)
}

// --------------------------------------

pub fn ecef2eci(gst :f64,x: f64, y: f64, z: f64) -> (f64,f64,f64){
    let arr = matmul3(transpose3(r3(gst)),[x,y,z]);
    (arr[0],arr[1],arr[2])

}

pub fn ecef2geodetic(x: f64, y: f64, z: f64) -> (f64,f64,f64){
    let major = wgs84().0;
    let minor = wgs84().1;

    let r = (x*x+y*y+z*z).sqrt();
    let e = (major*major-minor*minor).sqrt();
    let var = r*r-e*e;
    let u = (0.5*var+0.5*(var*var+4.0*e*e*z*z).sqrt()).sqrt();

    let q = (x*x+y*y).sqrt();
    let hu_e = (u*u+e*e).sqrt();
    let mut beta = (hu_e/u * z/q).atan();

    let eps =  ((minor * u - major * hu_e + e*e)*beta.sin())/
                    (major * hu_e/beta.cos() - e*e*beta.cos());
    beta += eps;

    let lat = (major/minor*beta.tan()).atan();
    let lon = y.atan2(x);

    let v1 = z - minor * beta.sin();
    let v2 = q - major*beta.cos();
    let alt;

    let inside = (x*x/major/major)+(y*y/major/major)+(z*z/minor/minor)<1.0;
    if inside {
        alt = -(v1*v1+v2*v2).sqrt();
    }else {
        alt = (v1*v1+v2*v2).sqrt();
    };

    (lat,lon,alt)
} 

pub fn ecef2enu(x: f64, y: f64, z: f64, lat0: f64, lon0: f64, alt0: f64)-> (f64,f64,f64){
    let (x0,y0,z0) = geodetic2ecef(lat0, lon0, alt0);
    let (e,n,u) = ecef2enuv(x-x0,y-y0,z-z0,lat0,lon0);
    (e,n,u)
}

pub fn ecef2enuv(u : f64, v: f64,w :f64, lat0: f64, lon0: f64) -> (f64,f64,f64){
    let t = lon0.cos() * u + lon0.sin() * v;
    let e = -lon0.sin() * u + lon0.cos() * v;
    let n = -lat0.sin() * t + lat0.cos() * w;
    let u = lat0.cos() * t + lat0.sin() * w;
    (e,n,u)
}

pub fn ecef2aer(x: f64, y: f64, z: f64, lat0: f64, lon0: f64, alt0: f64)-> (f64,f64,f64){
    let (e,n,u) = ecef2enu(x, y, z, lat0, lon0, alt0);
    let (az,el,slant_range) = enu2aer(e,n,u);

    (az,el,slant_range)
}

// --------------------------------------

pub fn eci2aer(gst :f64,x: f64, y: f64, z: f64, lat :f64, lon :f64, alt :f64)-> (f64,f64,f64){
    let (x,y,z) = eci2ecef(gst, x, y, z);
    let (az,el,slant_range) = ecef2aer(x, y, z, lat, lon, alt);
    (az,el,slant_range)
}

pub fn eci2ecef(gst :f64,x: f64, y: f64, z: f64)-> (f64,f64,f64){
    let arr = matmul3(r3(gst),[x,y,z]);
    (arr[0],arr[1],arr[2])
}

// --------------------------------------
pub fn matmul3(matrix: [f64;9],col:[f64;3])->[f64;3]{
    let out : [f64;3] = [ matrix[0]*col[0]+matrix[1]*col[1]+matrix[2]*col[2],
                        matrix[2]*col[0]+matrix[3]*col[1]+matrix[4]*col[2],
                        matrix[5]*col[0]+matrix[6]*col[1]+matrix[7]*col[2]];
    out
}
pub fn r3(x : f64) -> [f64; 9] {
    [x.cos(),x.sin(),0.0,-x.sin(),x.cos(),0.0,0.0,0.0,1.0]
}
pub fn transpose3(x: [f64;9]) -> [f64;9] {
    [x[0],x[3],x[6],x[1],x[4],x[7],x[2],x[5],x[8]]
}
pub fn wgs84() -> (f64,f64,f64,f64) {

    let major = 6378137.0;
    let flattening = 1.0/298.2572235630;
    let minor = major * (1.0 - flattening);
    let ecc_sq = ((major*major)-(minor*minor))/(major*major);
    
    (major,minor,flattening,ecc_sq)
}
pub fn get_radius_normal(lat: f64)->f64 {
    let (major,_,_,squared_eccentricity) = wgs84();
    major/((1.0-squared_eccentricity*lat.sin()*lat.sin()).sqrt())
}

pub fn deg2rad(x: f64) -> f64 {
    x/180.0*std::f64::consts::PI
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_geodetic() {

        let lat = deg2rad(30.14988205);
        let lon = deg2rad(91.38733072);
        let alt = 4031.0;

        let (x,y,z) = geodetic2ecef(lat,lon,alt);

        let xref = -1.337281037300386e+05;
        let yref = 5.521796910920261e+06;
        let zref = 3.186776473672415e+06;

        assert!((x-xref).abs()<1e-3);
        assert!((y-yref).abs()<1e-3);
        assert!((z-zref).abs()<1e-3);


    }


    
}


    
