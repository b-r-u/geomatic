use geomatic::laea;
use geomatic::{Point3035, Point4326};

fn main() {
    // Convert a WGS84 point to Point3035

    let p: Point4326 = Point4326::new(50.0, 5.0);
    let q: Point3035 = laea::forward_ref(p);

    println!("{} -> {}", p, q);
}
