// Nova Rust, a core library wrapping the Swiss Ephemeris library for astrological use.
// Copyright (C) 2024 Mike Verducci
/*This file is part of Nova Rust.

Nova Rust is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Nova Rust is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.
*/

pub mod file_utils;
pub mod spacetime_utils;
pub mod structs;
pub mod utils;

use crate::structs::{
    AngularityOrb, CoordinateSystemValues, Location, MeasuringFramework, Planet, SiderealContext,
};
use crate::utils::decimal_to_dms;
use anyhow::{anyhow, Result};
use chrono::{DateTime, Datelike, Days, NaiveDate, TimeDelta, TimeZone, Timelike, Utc};
use chrono_tz::Tz;
use libc::{c_char, c_int, c_longlong};
use phf::phf_map;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::Add;
use structs::{
    Angularities, CoordinateRange, HorizonForegroundLongitudes, NonHorizonForegroundLongitudes,
};

use reverse_geocoder::{ReverseGeocoder, SearchResult};
use std::ffi::{c_double, CString};
use std::{fmt, ptr};

use utils::{
    calculate_prime_vertical_longitude, calculate_right_ascension, clamp_results_to_range, format_two_decimals, get_opposite_geo_longitude, parse_angularity_longitude, parse_angularity_pvl, parse_angularity_ra
};

const SIDEREALMODE: c_int = 64 * 1024;
const CAMPANUS: c_int = 64;

const SUN_ORBITAL_HOURS: i32 = 8767;
const MOON_ORBITAL_HOURS: i32 = 656;

const ZODIAC: [&'static str; 12] = [
    "Ari", "Tau", "Gem", "Can", "Leo", "Vir", "Lib", "Sco", "Sag", "Cap", "Aqu", "Pis",
];

pub const PLANET_TO_INT: phf::Map<&'static str, i32> = phf_map! {
        "Sun" => 0,
        "Moon" => 1,
        "Mercury" => 2,
        "Venus" => 3,
        "Mars" => 4,
        "Jupiter" => 5,
        "Saturn" => 6,
        "Uranus" => 7,
        "Neptune" => 8,
        "Pluto" => 9,
        // minor planet number + 10,000 as SWE documentation instructs
        "Quaoar" => 60000,
        "Sedna" => 100377,
        "Orcus" => 100482,
        "Haumea" => 146108,
        "Eris" => 146199,
        "Makemake" => 146472,
        "Gonggong" => 235088,
};

pub const PLANET_NAMES: [&'static str; 17] = [
    "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto",
    "Quaoar", "Sedna", "Orcus", "Haumea", "Eris", "Makemake", "Gonggong",
];

fn body_is_tno(body_number: i32) -> bool {
    [60000, 100377, 100482, 146108, 146472, 146199, 235088].contains(&body_number)
}

fn body_is_malefic(body_number: i32) -> bool {
    [4, 6, 8].contains(&body_number)
}

fn body_is_benefic(body_number: i32) -> bool {
    [0, 3, 5, 7].contains(&body_number)
}

// Minor planet numbers (which don't have 10k added to them)
// Quaoar: 50000
// Orcus: 90482  (note - at around 950 km diameter, this is just below the 1k km cutoff)
// Haumea: 136108
// Makemake: 136472
// Gonggong: 225088

extern "C" {
    fn swe_set_ephe_path(path: *const c_char);
    fn swe_set_sid_mode(tjd_ut: f64, iflag: c_int, serr: *mut c_char) -> c_int;
    fn swe_close();
    fn swe_calc_ut(
        julian_day: c_double,
        body_number: c_int,
        flag: c_int,
        result: *mut c_double,
        err: *mut c_char,
    );
    fn swe_julday(
        year: c_int,
        month: c_int,
        day: c_int,
        hour: c_double,
        gregflag: c_int,
    ) -> c_double;
    fn swe_get_ayanamsa_ex_ut(
        julian_day: c_double,
        mode: c_int,
        svp: *mut c_double,
        err: *mut c_char,
    ) -> c_int;
    fn swe_houses_ex(
        julian_day: c_double,
        mode: c_int,
        latitude: c_double,
        longitude: c_double,
        house_system: c_int,
        first_cusp: *mut c_double,
        first_angle: *mut c_double,
    );
}

pub fn open_dll() {
    let current_dir = std::env::current_dir().unwrap();
    let path = current_dir
        .join("src")
        .join("swe")
        .join("ephemeris")
        .as_os_str()
        .to_string_lossy()
        .to_string();
    let c_str = CString::new(path).expect("CString::new failed");
    let c_str_ptr = c_str.as_ptr();
    unsafe {
        swe_set_ephe_path(c_str_ptr);
    }

    unsafe {
        swe_set_sid_mode(0.0, 0, ptr::null_mut());
    }
}

pub fn close_dll() {
    unsafe {
        swe_close();
    }
}

fn calculate_planets_ut(julian_day: f64, body_number: i32) -> Result<Vec<f64>> {
    unsafe {
        let mut swe_array: [f64; 6] = Default::default();
        let mut err_buffer: [c_char; 256] = [0; 256];

        swe_calc_ut(
            julian_day as c_double,
            body_number as c_int,
            (64 * 1024) as c_int,
            swe_array.as_mut_ptr(),
            err_buffer.as_mut_ptr(),
        );

        let float_array: Vec<f64> = swe_array.iter().cloned().collect();
        let err_string = error_string_from_buffer(&err_buffer);

        match err_string {
            Ok(e_string) => match e_string {
                Some(err) => Err(anyhow!(err)),
                None => Ok(float_array),
            },
            Err(e) => Err(anyhow!(e)),
        }
    }
}

fn get_orb(coord1: f64, coord2: f64, low_bound: f64, high_bound: f64) -> Option<f64> {
    let aspect_average = (low_bound + high_bound) / 2.0;
    let aspect = f64::abs(coord1 - coord2);

    // If one longitude is near 360º and the other is near 0º
    let aspect360 = f64::abs(aspect - 360.0);

    if low_bound <= aspect && aspect <= high_bound {
        return match low_bound as i32 {
            0 => Some(aspect),
            _ => Some(f64::abs(aspect - aspect_average)),
        };
    }

    if low_bound <= aspect360 && aspect360 <= high_bound {
        return match low_bound as i32 {
            0 => Some(aspect360),
            _ => Some(f64::abs(aspect360 - aspect_average)),
        };
    }

    None
}

fn get_orb_signed(
    base_coord: f64,
    reference_coord: f64,
    low_bound: f64,
    high_bound: f64,
) -> Option<f64> {
    let aspect_average = (low_bound + high_bound) / 2.0;
    let raw_orb = base_coord - reference_coord;
    let aspect = f64::abs(raw_orb);

    let aspect_average = (low_bound + high_bound) / 2.0;
    let aspect = f64::abs(base_coord - reference_coord);

    // If one longitude is near 360º and the other is near 0º
    let aspect360 = f64::abs(aspect - 360.0);

    if low_bound <= aspect && aspect <= high_bound {
        return match low_bound as i32 {
            0 => match raw_orb < 0.0 {
                true => Some(aspect * -1.0),
                false => Some(aspect),
            },
            _ => match raw_orb < 0.0 {
                true => Some(f64::abs(aspect - aspect_average) * -1.0),
                false => Some(f64::abs(aspect - aspect_average)),
            },
        };
    }

    if low_bound <= aspect360 && aspect360 <= high_bound {
        return match low_bound as i32 {
            0 => match raw_orb < 0.0 {
                true => Some(aspect360 * -1.0),
                false => Some(aspect360),
            },
            _ => match raw_orb < 0.0 {
                true => Some(f64::abs(aspect360 - aspect_average) * -1.0),
                false => Some(f64::abs(aspect360 - aspect_average)),
            },
        };
    }

    None
}

fn calculate_cusps_and_angles(
    julian_day: f64,
    longitude: f64,
    latitude: f64,
) -> ([f64; 13], [f64; 10]) {
    let mut cusps: [f64; 13] = Default::default();
    let mut angles: [f64; 10] = Default::default();

    unsafe {
        swe_houses_ex(
            julian_day,
            SIDEREALMODE,
            latitude,
            longitude,
            CAMPANUS,
            cusps.as_mut_ptr(),
            angles.as_mut_ptr(),
        );
    }

    // cusps are numbered
    // angles are: asc, mc, armc, vertex, "equatorial ascendant," and stuff we don't care about
    (cusps, angles)
}

fn get_julian_day(utc_dt: DateTime<Utc>) -> f64 {
    let decimal_hour = convert_dms_to_decimal(
        utc_dt.hour() as f64,
        utc_dt.minute() as f64,
        utc_dt.second() as f64,
    );

    unsafe {
        let julian_day = swe_julday(
            utc_dt.year() as c_int,
            utc_dt.month() as c_int,
            utc_dt.day() as c_int,
            decimal_hour,
            1,
        );

        julian_day
    }
}

fn convert_dms_to_decimal(degrees: f64, minutes: f64, seconds: f64) -> c_double {
    degrees + minutes / 60.0 + seconds / 3600.0
}

fn error_string_from_buffer(err_buffer: &[c_char; 256]) -> anyhow::Result<Option<String>> {
    let err_cstring: Vec<u8> = err_buffer.iter().map(|&c| c as u8).collect();
    let err_string = String::from_utf8(err_cstring)?;
    let trimmed_err_string = err_string.trim_matches(char::from(0)).to_string();
    match trimmed_err_string.len() {
        0 => Ok(None),
        _ => Ok(Some(trimmed_err_string)),
    }
}

fn get_julian_day_0_gmt(utc_dt: chrono::DateTime<chrono::Utc>) -> f64 {
    unsafe {
        let julian_day = swe_julday(
            utc_dt.year() as c_int,
            utc_dt.month() as c_int,
            utc_dt.day() as c_int,
            0.0,
            1,
        );

        julian_day
    }
}

fn calculate_lst(dt: DateTime<Utc>, longitude: f64) -> f64 {
    let julian_day_0_gmt = get_julian_day_0_gmt(dt);

    let universal_time =
        convert_dms_to_decimal(dt.hour() as f64, dt.minute() as f64, dt.second() as f64);

    let sidereal_time_at_midnight_julian_day = (julian_day_0_gmt - 2451545.0) / 36525.0;

    let greenwich_sidereal_time = 6.697374558
        + 2400.051336 * sidereal_time_at_midnight_julian_day
        + 0.000024862 * sidereal_time_at_midnight_julian_day.powi(2)
        + universal_time * 1.0027379093;

    let local_sidereal_time = (greenwich_sidereal_time + (longitude / 15.0)).rem_euclid(24.0);

    if local_sidereal_time < 0.0 {
        local_sidereal_time + 24.0
    } else {
        local_sidereal_time
    }
}

fn get_obliquity(julian_day: f64) -> Result<f64> {
    let obliquity_array = calculate_planets_ut(julian_day, -1)?;
    Ok(obliquity_array[0])
}

pub fn create_sidereal_context(dt: DateTime<Tz>, location: Location) -> Result<SiderealContext> {
    let utc_dt = dt.with_timezone(&Utc);
    let julian_day = get_julian_day(utc_dt);

    let svp = get_svp(julian_day)?;
    let lst = calculate_lst(utc_dt, location.longitude);
    let ramc = lst * 15.0;
    let obliquity = get_obliquity(julian_day)?;

    Ok(SiderealContext {
        lst,
        ramc,
        obliquity,
        svp,
        datetime: utc_dt,
        julian_day,
        location,
    })
}

fn get_svp(julian_day: f64) -> Result<f64> {
    let mut swe_array: [f64; 1] = Default::default();
    let mut err_buffer: [c_char; 256] = [0; 256];

    unsafe {
        let is_err = swe_get_ayanamsa_ex_ut(
            julian_day,
            SIDEREALMODE,
            swe_array.as_mut_ptr(),
            err_buffer.as_mut_ptr(),
        );
        match is_err < 0 {
            true => Err(anyhow!("Could not calculate SVP")),
            false => Ok(30.0 - swe_array[0]),
        }
    }
}

// For each degree of latitude between 24 and 50
// for each degree of longitude between -66 and -125
// (total of 1,534 locations)
// get LST for that location (really, precess the natal planets to there) and determine angular natal planets.
// If there are any, geocode that lat/long and add the place name to a list.
// at the end, print.

// If we do this globally: extend latitude to from -66 to 66
// and extend longitude from -170 to 180
// This brings us up to 43,050 locations.

pub fn find_next_solunar_utc_dt(
    body_number: i32,
    radix_position: f64,
    harmonic: i8,
    time_reverse: bool,
    base_dt: DateTime<Utc>,
) -> anyhow::Result<DateTime<Utc>> {
    let mut closest_approach_orb = 359.0;

    let mut test_dt = base_dt
        .checked_add_days(Days::new(1))
        .ok_or_else(|| anyhow::Error::msg("No date found"))?
        .with_hour(0)
        .ok_or_else(|| anyhow::Error::msg("Cannot add hour 0 to date"))?;
    let mut julian_day = get_julian_day(test_dt);
    let mut transiting_data = calculate_planets_ut(julian_day, body_number)?;

    // Search by day - will eventually do as a binary search
    loop {
        match get_orb_signed(radix_position, transiting_data[0], 0.0, 5.0) {
            Some(current_orb) => {
                // If we were closer before and went past
                if current_orb < 0.0 {
                    test_dt = test_dt
                        .checked_sub_days(Days::new(1))
                        .ok_or_else(|| anyhow::Error::msg("Could not subtract days"))?;
                    julian_day = get_julian_day(test_dt);
                    transiting_data = calculate_planets_ut(julian_day, body_number)?;
                    break;
                }
                test_dt = test_dt
                    .checked_add_days(Days::new(1))
                    .ok_or_else(|| anyhow::Error::msg("Could not add days"))?;
                julian_day = get_julian_day(test_dt);
                transiting_data = calculate_planets_ut(julian_day, body_number)?;
            }
            None => {
                test_dt = test_dt
                    .checked_add_days(Days::new(1))
                    .ok_or_else(|| anyhow::Error::msg("Could not add days"))?;
                julian_day = get_julian_day(test_dt);
                transiting_data = calculate_planets_ut(julian_day, body_number)?;
            }
        };
    }

    // Search by second - will eventually do as a binary search
    loop {
        match get_orb_signed(radix_position, transiting_data[0], 0.0, 5.0) {
            Some(current_orb) => {
                // If we were closer before and went past
                if current_orb < 0.0 {
                    test_dt = test_dt.add(
                        TimeDelta::new(-1, 0)
                            .ok_or_else(|| anyhow::Error::msg("Could not add timedelta"))?,
                    );
                    break;
                }
                test_dt = test_dt.add(
                    TimeDelta::new(1, 0)
                        .ok_or_else(|| anyhow::Error::msg("Could not add timedelta"))?,
                );
                julian_day = get_julian_day(test_dt);
                transiting_data = calculate_planets_ut(julian_day, body_number)?;
            }
            None => {
                test_dt = test_dt.add(
                    TimeDelta::new(1, 0)
                        .ok_or_else(|| anyhow::Error::msg("Could not add timedelta"))?,
                );
                julian_day = get_julian_day(test_dt);
                transiting_data = calculate_planets_ut(julian_day, body_number)?;
            }
        };
    }

    Ok(test_dt)
}

fn sorted_angular_planets_by_orb(angular_planets: Vec<AngularityOrb>) -> Angularities {
    let mut sorted_angular_planets = angular_planets.clone();
    sorted_angular_planets.sort_by(|a, b| {
        let a_orb = a.orb.0 as f64 + a.orb.1 as f64 / 60.0;
        let b_orb = b.orb.0 as f64 + b.orb.1 as f64 / 60.0;
        a_orb.partial_cmp(&b_orb).unwrap()
    });

    // Take only the closest orb for each planet
    let planets =
        sorted_angular_planets
            .iter()
            .fold(Vec::<AngularityOrb>::new(), |mut acc, angularity| {
                if !acc.iter().any(|a| a.planet == angularity.planet) {
                    let x = AngularityOrb {
                        planet: angularity.planet.clone(),
                        framework: angularity.framework.clone(),
                        orb: angularity.orb,
                    };
                    acc.push(x);
                }
                acc
            });

    Angularities { planets }
}

fn sort_list_of_angularities_by_closest_orb(
    angularities: HashMap<String, Angularities>,
) -> Vec<(String, Angularities)> {
    let mut sorted_angularities: Vec<(String, Angularities)> = angularities.into_iter().collect();
    sorted_angularities.sort_by(|a, b| {
        let a_orb = a.1.planets[0].orb.0 as f64 + a.1.planets[0].orb.1 as f64 / 60.0;
        let b_orb = b.1.planets[0].orb.0 as f64 + b.1.planets[0].orb.1 as f64 / 60.0;
        a_orb.partial_cmp(&b_orb).unwrap()
    });

    sorted_angularities
}

pub fn find_longitudes_for_planet_on_meridian_axis(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    precision: f64,
) -> anyhow::Result<Vec<f64>> {
    let julian_day_utc = get_julian_day(utc_dt);
    let geo_latitude = 0.0;

    let mut orb = 0.0;

    // find MC; IC is just on the opposite longitude
    let mut upper = 180.0;
    let mut lower = -180.0;
    let mut geo_longitude: f64 = 0.0;
    loop {
        if upper < lower {
            break;
        }

        geo_longitude = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude) * 15.0;
        let obliquity = get_obliquity(julian_day_utc)?;
        let svp = get_svp(julian_day_utc)?;

        let pvl = calculate_prime_vertical_longitude(
            planet_coords.longitude,
            planet_coords.latitude,
            ramc,
            obliquity,
            svp,
            geo_latitude,
        );

        orb = get_orb_signed(pvl, 270.0, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude - 0.01,
            false => lower = geo_longitude + 0.01,
        };
    }

    Ok(clamp_results_to_range(
        Vec::from([geo_longitude, get_opposite_geo_longitude(geo_longitude)]),
        range,
    ))
}

pub fn find_longitudes_for_planet_square_meridian(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    precision: f64,
) -> anyhow::Result<Vec<f64>> {
    let julian_day_utc = get_julian_day(utc_dt);
    let geo_latitude = 0.0;

    let mut orb = 0.0;

    // find EP; WP is just on the opposite longitude
    let mut upper = 180.0;
    let mut lower = -180.0;
    let mut geo_longitude: f64 = 0.0;
    loop {
        if upper < lower {
            break;
        }

        geo_longitude = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude) * 15.0;
        let obliquity = get_obliquity(julian_day_utc)?;
        let svp = get_svp(julian_day_utc)?;

        let (_, angles) = calculate_cusps_and_angles(julian_day_utc, geo_longitude, geo_latitude);
        let mc = angles[1];

        let rising_square_mc = (mc - 90.0).rem_euclid(360.0);

        orb = get_orb_signed(planet_coords.longitude, rising_square_mc, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude - 0.01,
            false => lower = geo_longitude + 0.01,
        };
    }

    Ok(clamp_results_to_range(
        Vec::from([geo_longitude, get_opposite_geo_longitude(geo_longitude)]),
        range,
    ))
}

pub fn find_longitudes_for_planet_square_ramc(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    precision: f64,
) -> anyhow::Result<Vec<f64>> {
    let julian_day_utc = get_julian_day(utc_dt);
    let geo_latitude = 0.0;

    let mut orb = 0.0;

    let obliquity = get_obliquity(julian_day_utc)?;
    let svp = get_svp(julian_day_utc)?;
    let ra = calculate_right_ascension(
        planet_coords.longitude,
        planet_coords.latitude,
        svp,
        obliquity,
    );

    // find one of EP-a or WP-a; the other is just on the opposite longitude
    let mut upper = 180.0;
    let mut lower = -180.0;
    let mut geo_longitude: f64 = 0.0;
    loop {
        if upper < lower {
            break;
        }

        geo_longitude = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude) * 15.0;

        let ep_a = (ramc + 90.0).rem_euclid(360.0);

        orb = get_orb_signed(ra, ep_a, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude - 0.01,
            false => lower = geo_longitude + 0.01,
        };
    }

    Ok(clamp_results_to_range(
        Vec::from([geo_longitude, get_opposite_geo_longitude(geo_longitude)]),
        range,
    ))
}

pub fn find_longitudes_for_planet_on_horizon(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    latitude: i32,
    precision: f64,
) -> anyhow::Result<Vec<f64>> {
    let julian_day_utc = get_julian_day(utc_dt);

    let obliquity = get_obliquity(julian_day_utc)?;
    let svp = get_svp(julian_day_utc)?;

    let mut upper = 180.0;
    let mut lower = -180.0;
    let mut geo_longitude_1: f64 = 0.0;
    let mut orb = 0.0;

    let mut results: Vec<f64> = Vec::new();

    // Asc or Dsc
    loop {
        if upper < lower {
            break;
        }

        geo_longitude_1 = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude_1) * 15.0;

        let pvl = calculate_prime_vertical_longitude(
            planet_coords.longitude,
            planet_coords.latitude,
            ramc,
            obliquity,
            svp,
            latitude as f64,
        );

        orb = get_orb_signed(pvl, 0.0, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude_1 - 0.01,
            false => lower = geo_longitude_1 + 0.01,
        };
    }

    if geo_longitude_1 > range.min_longitude as f64 && geo_longitude_1 < range.max_longitude as f64
    {
        results.push(geo_longitude_1);
    }

    upper = 180.0;
    lower = -180.0;

    // Exclude the first longitude so we don't accidentally find it again
    match geo_longitude_1 > 0.0 {
        true => upper = geo_longitude_1 - 1.0,
        false => lower = geo_longitude_1 + 1.0,
    };

    let mut geo_longitude_2 = 0.0;
    loop {
        if upper < lower {
            break;
        }

        geo_longitude_2 = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude_2) * 15.0;

        let pvl = calculate_prime_vertical_longitude(
            planet_coords.longitude,
            planet_coords.latitude,
            ramc,
            obliquity,
            svp,
            latitude as f64,
        );

        orb = get_orb_signed(pvl, 180.0, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude_2 - 0.01,
            false => lower = geo_longitude_2 + 0.01,
        };
    }
    // Find the other of Asc or Dsc

    if geo_longitude_2 > range.min_longitude as f64 && geo_longitude_2 < range.max_longitude as f64
    {
        results.push(geo_longitude_2);
    }

    Ok(results)
}

pub fn find_horizon_conjunctions_within_latitude_range(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    latitude_band_size: i32,
    precision: f64,
) -> anyhow::Result<Vec<(i32, f64)>> {
    let julian_day_utc = get_julian_day(utc_dt);

    // latitude, asc or dsc
    let mut conjunctions: Vec<(i32, f64)> = Vec::new();
    for latitude in (range.min_latitude..=range.max_latitude).step_by(latitude_band_size as usize) {
        find_longitudes_for_planet_on_horizon(utc_dt, planet_coords, range, latitude, precision)?
            .into_iter()
            .for_each(|longitude| {
                conjunctions.push((latitude, longitude));
            });
    }

    Ok(conjunctions)
}

pub fn find_longitudes_for_planet_square_horizon(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    latitude: i32,
    precision: f64,
) -> anyhow::Result<Vec<f64>> {
    let julian_day_utc = get_julian_day(utc_dt);

    let mut orb = 0.0;

    let mut results: Vec<f64> = Vec::new();

    // find one of Zenith or Nadir
    let mut upper = 180.0;
    let mut lower = -180.0;
    let mut geo_longitude_1: f64 = 0.0;
    loop {
        if upper < lower {
            break;
        }

        geo_longitude_1 = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude_1) * 15.0;
        let obliquity = get_obliquity(julian_day_utc)?;
        let svp = get_svp(julian_day_utc)?;

        let (_, angles) =
            calculate_cusps_and_angles(julian_day_utc, geo_longitude_1, latitude as f64);
        let asc = angles[0];

        let culminating_square = (asc + 90.0).rem_euclid(360.0);

        orb = get_orb_signed(planet_coords.longitude, culminating_square, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude_1 - 0.01,
            false => lower = geo_longitude_1 + 0.01,
        };
    }

    if geo_longitude_1 > range.min_longitude as f64 && geo_longitude_1 < range.max_longitude as f64
    {
        results.push(geo_longitude_1);
    }

    // find the other of of Zenith or Nadir
    upper = 180.0;
    lower = -180.0;

    match geo_longitude_1 > 0.0 {
        true => upper = geo_longitude_1 - 1.0,
        false => lower = geo_longitude_1 + 1.0,
    };

    let mut geo_longitude_2 = 0.0;
    loop {
        if upper < lower {
            break;
        }

        geo_longitude_2 = (upper + lower) / 2.0;

        let ramc = calculate_lst(utc_dt, geo_longitude_2) * 15.0;
        let obliquity = get_obliquity(julian_day_utc)?;
        let svp = get_svp(julian_day_utc)?;

        let (_, angles) =
            calculate_cusps_and_angles(julian_day_utc, geo_longitude_2, latitude as f64);
        let asc = angles[0];

        let anti_culminating_square = (asc - 90.0).rem_euclid(360.0);

        orb = get_orb_signed(planet_coords.longitude, anti_culminating_square, 0.0, 360.0).unwrap();

        // Normalize to 180º
        orb = match orb as i32 {
            -360..=-180 => orb + 360.0,
            180..=360 => orb - 360.0,
            _ => orb,
        };

        if f64::abs(orb) < precision {
            break;
        };

        match orb < 0.0 {
            true => upper = geo_longitude_2 - 0.01,
            false => lower = geo_longitude_2 + 0.01,
        };
    }

    if geo_longitude_2 > range.min_longitude as f64 && geo_longitude_2 < range.max_longitude as f64
    {
        results.push(geo_longitude_2);
    }

    Ok(results)
}

/// Returns latitude, and the two longitudes for the square
pub fn find_horizon_squares_within_latitude_range(
    utc_dt: DateTime<Utc>,
    planet_coords: &CoordinateSystemValues,
    range: &CoordinateRange,
    latitude_band_size: i32,
    precision: f64,
) -> anyhow::Result<Vec<(i32, f64)>> {
    let julian_day_utc = get_julian_day(utc_dt);

    // latitude, asc, dsc
    let mut conjunctions: Vec<(i32, f64)> = Vec::new();
    for latitude in (range.min_latitude..=range.max_latitude).step_by(latitude_band_size as usize) {
        find_longitudes_for_planet_square_horizon(
            utc_dt,
            planet_coords,
            range,
            latitude,
            precision,
        )?
        .into_iter()
        .for_each(|longitude| {
            conjunctions.push((latitude, longitude));
        });
    }

    Ok(conjunctions)
}

pub fn foreground_coordinates_in_solunar_return(
    radix_dt: DateTime<Tz>,
    target_dt: DateTime<Tz>,
    range: crate::structs::CoordinateRange,
    allowed_body_numbers: Vec<i32>,
    body_number: i32,
    harmonic: i32,
) -> anyhow::Result<(NonHorizonForegroundLongitudes, HorizonForegroundLongitudes)> {
    let mut angular_locations: HashMap<String, Vec<AngularityOrb>> = HashMap::new();

    let radix_utc = radix_dt.with_timezone(&Utc);
    let radix_julian_day = get_julian_day(radix_utc);

    let radix_positions = calculate_planets_ut(radix_julian_day, body_number)?;
    let target_utc_dt = target_dt.with_timezone(&Utc);

    let solunar_dt = find_next_solunar_utc_dt(
        body_number,
        radix_positions[0],
        harmonic as i8,
        false,
        target_utc_dt,
    )?;

    let solunar_julian_day = get_julian_day(solunar_dt);

    let svp = get_svp(solunar_julian_day)?;
    let obliquity = get_obliquity(solunar_julian_day)?;

    let radix_positions: anyhow::Result<Vec<Planet>> = PLANET_TO_INT
        .entries()
        .into_iter()
        .map(
            |(name, &body_number)| match calculate_planets_ut(radix_julian_day, body_number) {
                Ok(values) => Ok(Planet {
                    name: name.to_string(),
                    body_number,
                    coordinates: {
                        CoordinateSystemValues {
                            longitude: values[0],
                            latitude: values[1],
                            speed_in_longitude: values[2],
                            right_ascension: 0.0,
                            prime_vertical_longitude: 0.0,
                            azimuth: 0.0,
                            meridian_longitude: 0.0,
                        }
                    },
                }),
                Err(e) => Err(e),
            },
        )
        .collect();

    if radix_positions.is_err() {
        return Err(radix_positions.err().unwrap());
    }

    let solunar_positions: anyhow::Result<Vec<Planet>> = PLANET_TO_INT
        .entries()
        .into_iter()
        .map(
            |(name, &body_number)| match calculate_planets_ut(solunar_julian_day, body_number) {
                Ok(values) => Ok(Planet {
                    name: name.to_string(),
                    body_number,
                    coordinates: {
                        CoordinateSystemValues {
                            longitude: values[0],
                            latitude: values[1],
                            speed_in_longitude: values[2],
                            right_ascension: 0.0,
                            prime_vertical_longitude: 0.0,
                            azimuth: 0.0,
                            meridian_longitude: 0.0,
                        }
                    },
                }),
                Err(e) => Err(e),
            },
        )
        .collect();

    // latitude: (planet, longitude)
    let mut latitude_bands: HashMap<i32, Vec<(String, f64)>> = HashMap::new();
    let mut longitude_lines: HashMap<String, Vec<f64>> = HashMap::new();

    for planet in radix_positions.as_ref().unwrap() {
        if !allowed_body_numbers.contains(&planet.body_number) {
            continue;
        }

        let mut planet_name_with_context = String::from("(R) ");
        planet_name_with_context.push_str(planet.name.as_str());

        longitude_lines.insert(planet_name_with_context.clone(), Vec::new());

        let meridian_conj_longitudes = find_longitudes_for_planet_on_meridian_axis(
            solunar_dt,
            &planet.coordinates,
            &range,
            0.01,
        )?;

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .extend(meridian_conj_longitudes);

        let meridian_square_longitudes = find_longitudes_for_planet_square_meridian(
            solunar_dt,
            &planet.coordinates,
            &range,
            0.01,
        )?;

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .extend(meridian_square_longitudes);

        let ramc_square_longitudes =
            find_longitudes_for_planet_square_ramc(solunar_dt, &planet.coordinates, &range, 0.01)?;

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .extend(ramc_square_longitudes);

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .sort_by(|a, b| a.partial_cmp(b).unwrap());

        let horizon_conj_points = find_horizon_conjunctions_within_latitude_range(
            solunar_dt,
            &planet.coordinates,
            &range,
            5,
            0.01,
        )?;

        horizon_conj_points.into_iter().for_each(|(latitude, longitude)| {
            latitude_bands.entry(latitude).or_insert(Vec::new());
            let mut planet_name_with_context = String::from("(R) ");
            planet_name_with_context.push_str(planet.name.as_str());

            latitude_bands
                .get_mut(&latitude)
                .unwrap()
                .push((planet_name_with_context.clone(), longitude));
        });

        let horizon_square_points = find_horizon_squares_within_latitude_range(
            solunar_dt,
            &planet.coordinates,
            &range,
            5,
            0.01,
        )?;

        horizon_square_points.into_iter().for_each(|(latitude, longitude)| {
            latitude_bands.entry(latitude).or_insert(Vec::new());
            let mut planet_name_with_context = String::from("(R) ");
            planet_name_with_context.push_str(planet.name.as_str());

            latitude_bands
                .get_mut(&latitude)
                .unwrap()
                .push((planet_name_with_context.clone(), longitude));
        });
    }

    for planet in solunar_positions.as_ref().unwrap() {
        if !allowed_body_numbers.contains(&planet.body_number) {
            continue;
        }

        let mut planet_name_with_context = String::from("(T) ");
        planet_name_with_context.push_str(planet.name.as_str());

        longitude_lines.insert(planet_name_with_context.clone(), Vec::new());

        let meridian_conj_longitudes = find_longitudes_for_planet_on_meridian_axis(
            solunar_dt,
            &planet.coordinates,
            &range,
            0.01,
        )?;

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .extend(meridian_conj_longitudes);

        let meridian_square_longitudes = find_longitudes_for_planet_square_meridian(
            solunar_dt,
            &planet.coordinates,
            &range,
            0.01,
        )?;

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .extend(meridian_square_longitudes);

        let ramc_square_longitudes =
            find_longitudes_for_planet_square_ramc(solunar_dt, &planet.coordinates, &range, 0.01)?;

        longitude_lines
            .get_mut(&planet_name_with_context)
            .unwrap()
            .extend(ramc_square_longitudes);

        let horizon_conj_points = find_horizon_conjunctions_within_latitude_range(
            solunar_dt,
            &planet.coordinates,
            &range,
            5,
            0.01,
        )?;

        horizon_conj_points.into_iter().for_each(|(latitude, longitude)| {
            latitude_bands.entry(latitude).or_insert(Vec::new());
            let mut planet_name_with_context = String::from("(T) ");
            planet_name_with_context.push_str(planet.name.as_str());

            latitude_bands
                .get_mut(&latitude)
                .unwrap()
                .push((planet_name_with_context.clone(), longitude));
        });

        let horizon_square_points = find_horizon_squares_within_latitude_range(
            solunar_dt,
            &planet.coordinates,
            &range,
            5,
            0.01,
        )?;

        horizon_square_points.into_iter().for_each(|(latitude, longitude)| {
            latitude_bands.entry(latitude).or_insert(Vec::new());
            let mut planet_name_with_context = String::from("(T) ");
            planet_name_with_context.push_str(planet.name.as_str());

            latitude_bands
                .get_mut(&latitude)
                .unwrap()
                .push((planet_name_with_context.clone(), longitude));
        });
    }

    let mut longitude_vec = NonHorizonForegroundLongitudes {
        planets: Vec::new(),
    };

    for (planet_identifier, longitudes) in longitude_lines.iter() {
        let mut vec = longitudes.clone();
        vec.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let stringified_vec = vec.into_iter().map(|f| format_two_decimals(f)).collect::<Vec<String>>();
        longitude_vec.planets.push((planet_identifier.clone(), stringified_vec));
    }

    let mut latitude_vec = HorizonForegroundLongitudes {
        planets: Vec::new(),
    };

    for (&latitude, value) in latitude_bands.iter() {
        let mut vec = value.clone();
        vec.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let stringified_vec = vec.into_iter().map(|(key, f)| (key, format_two_decimals(f))).collect::<Vec<(String, String)>>();
        latitude_vec.planets.push((latitude, stringified_vec));
    }

    latitude_vec.planets.sort_unstable_by_key(|v| v.0);

    Ok((longitude_vec, latitude_vec))
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn t() {
        open_dll();
        let dt =
            spacetime_utils::string_to_naive_datetime("12/20/1989 10:20pm".to_string()).unwrap();
        let radix_dt = spacetime_utils::naive_to_local_tz(dt, -74.09184, 40.90345).unwrap();
        let target_dt_naive =
            spacetime_utils::string_to_naive_datetime("11/01/2024 10:00am".to_string()).unwrap();
        let target_dt =
            spacetime_utils::naive_to_local_tz(target_dt_naive, -74.09184, 40.90345).unwrap();

        let (longitudes, bands) = foreground_coordinates_in_solunar_return(
            radix_dt,
            target_dt,
            structs::CoordinateRange {
                min_longitude: -120,
                max_longitude: 120,
                min_latitude: 29,
                max_latitude: 49,
            },
            vec![3, 5, 7],
            0,
            1,
        )
        .unwrap();

        println!("{:?}", longitudes);
        println!("");
        println!("");
        println!("{:?}", bands);

        close_dll();
    }
}
