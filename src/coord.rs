use std::fmt;
use std::marker::PhantomData;
use crate::laea::LaeaParams;
use crate::ellipsoid::{Ellipsoid, GRS1980};

/// Coordinate reference system
pub trait Crs {
    const SRID: Option<u32>;
}
pub trait CrsLaea: Crs {
    type Ellipsoid: Ellipsoid;
    const PARAMS: LaeaParams;
}

#[derive(Copy, Clone, Debug)]
pub enum Crs4326 {}

impl Crs for Crs4326 {
    const SRID: Option<u32> = Some(4326);
}

/*
   PROJCS["ETRS89 / LAEA Europe",
    GEOGCS["ETRS89",
        DATUM["European_Terrestrial_Reference_System_1989",
            SPHEROID["GRS 1980",6378137,298.257222101,
                AUTHORITY["EPSG","7019"]],
            TOWGS84[0,0,0,0,0,0,0],
            AUTHORITY["EPSG","6258"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.0174532925199433,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4258"]],
    PROJECTION["Lambert_Azimuthal_Equal_Area"],
    PARAMETER["latitude_of_center",52],
    PARAMETER["longitude_of_center",10],
    PARAMETER["false_easting",4321000],
    PARAMETER["false_northing",3210000],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AUTHORITY["EPSG","3035"]]
*/

#[derive(Copy, Clone, Debug)]
pub enum Crs3035 {}

impl Crs for Crs3035 {
    const SRID: Option<u32> = Some(3035);
}

impl CrsLaea for Crs3035 {
    type Ellipsoid = GRS1980;

    const PARAMS: LaeaParams =
        LaeaParams {
            center_lat: 52.0,
            center_lon: 10.0,
            false_easting: 4321000.0,
            false_northing: 3210000.0,
        };
}

#[derive(Copy, Clone, Debug)]
pub struct Point<C: Crs> {
    pub coords: (f64, f64),
    _phantom: PhantomData<C>,
}

impl<C: Crs> Point<C> {
    pub fn new(a: f64, b: f64) -> Self {
        Self { coords: (a, b), _phantom: PhantomData }
    }
}

impl<C: Crs> fmt::Display for Point<C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match C::SRID {
            Some(srid) => write!(f, "Point({}, {}, srid={})", self.coords.0, self.coords.1, srid),
            None => write!(f, "Point({}, {}, srid=unknown)", self.coords.0, self.coords.1),
        }
    }
}

pub type Point4326 = Point<Crs4326>;

impl Point4326 {
    pub fn lat(&self) -> f64 {
        self.coords.0
    }

    pub fn lon(&self) -> f64 {
        self.coords.1
    }
}

pub type Point3035 = Point<Crs3035>;

impl Point3035 {
    pub fn east(&self) -> f64 {
        self.coords.0
    }

    pub fn north(&self) -> f64 {
        self.coords.1
    }
}
