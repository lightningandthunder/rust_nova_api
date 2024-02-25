use std::{fmt, path::Display};

use chrono::{
    DateTime, Datelike, Days, LocalResult, NaiveDate, NaiveDateTime, TimeDelta, TimeZone, Timelike,
    Utc,
};
use chrono_tz::Tz;

#[derive(Debug, PartialEq)]
pub struct CoordinateSystemValues {
    pub longitude: f64,
    pub latitude: f64,
    pub speed_in_longitude: f64,
    pub right_ascension: f64,
    pub prime_vertical_longitude: f64,
    pub azimuth: f64,
    pub meridian_longitude: f64, // This is experimental
}

#[derive(Debug, PartialEq)]
pub struct Planet {
    pub name: String,
    pub body_number: i32,
    pub coordinates: CoordinateSystemValues,
}

#[repr(C)]
pub union FlexibleDatetime {
    pub utc_datetime: DateTime<Utc>,
    pub local_datetime: DateTime<Tz>,
}

impl std::fmt::Debug for FlexibleDatetime {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        unsafe {
            match self {
                FlexibleDatetime { utc_datetime } => write!(f, "{}", utc_datetime),
                FlexibleDatetime { local_datetime } => write!(f, "{}", local_datetime),
            }
        }
    }
}

#[derive(Debug)]
pub struct Location {
    pub longitude: f64,
    pub latitude: f64,
    pub datetime: FlexibleDatetime,
    pub name: String,
}

#[derive(Debug, Clone)]
pub enum MeasuringFramework {
    Longitude,
    RightAscension,
    PrimeVerticalLongitude,
    Azimuth,
    MeridianLongitude,
}

impl fmt::Display for MeasuringFramework {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

#[derive(Debug, Clone)]
pub struct AngularityOrb {
    pub planet: String,
    pub framework: MeasuringFramework,
    pub orb: (i32, i32),
}

impl AngularityOrb {
    pub fn pretty_print(&self) -> String {
        format!(
            "{} {}° {}' ({})",
            self.planet, self.orb.0, self.orb.1, self.framework
        )
    }
}

impl fmt::Display for AngularityOrb {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} {}° {}' ({})",
            self.planet, self.orb.0, self.orb.1, self.framework,
        )
    }
}

pub struct Angularities {
    pub planets: Vec<AngularityOrb>,
}

impl fmt::Display for Angularities {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let planet_string = self
            .planets
            .iter()
            .map(|p| p.pretty_print())
            .collect::<Vec<String>>()
            .join(", ");
        write!(f, "{}", planet_string)
    }
}

pub struct AngularLocation {
    location: Location,
    planets: Vec<AngularityOrb>,
}

impl fmt::Display for AngularLocation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let planet_string = self
            .planets
            .iter()
            .map(|p| p.pretty_print())
            .collect::<Vec<String>>()
            .join(", ");
        write!(f, "{}: {}", self.location.name, planet_string)
    }
}

#[derive(Debug)]
pub struct SiderealContext {
    pub lst: f64,
    pub ramc: f64,
    pub obliquity: f64,
    pub svp: f64,
    pub datetime: DateTime<Utc>,
    pub julian_day: f64,
    pub location: Location,
}

pub struct CoordinateRange {
    pub min_latitude: i32,
    pub max_latitude: i32,
    pub min_longitude: i32,
    pub max_longitude: i32,
}
