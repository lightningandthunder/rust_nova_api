mod utils;

use chrono::{NaiveDate, TimeZone};
use nova_rust::{
    angular_precessed_planets_in_range, close_dll, create_sidereal_context, open_dll,
    CoordinateRange, FlexibleDatetime, Location,
};

fn main() {
    open_dll();
    let local = NaiveDate::from_ymd_opt(1989, 12, 20)
        .unwrap()
        .and_hms_opt(22, 20, 0)
        .unwrap();
    let radix_dt = chrono_tz::EST.from_local_datetime(&local).unwrap();

    let loc = Location {
        longitude: -74.1166,
        latitude: 40.97,
        datetime: FlexibleDatetime {
            local_datetime: radix_dt,
        },
        name: String::from("mv"),
    };

    let target = NaiveDate::from_ymd_opt(2023, 12, 21)
        .unwrap()
        .and_hms_opt(15, 29, 0)
        .unwrap();
    let target_dt = chrono_tz::EST.from_local_datetime(&target).unwrap();

    // For each degree of latitude between 24 and 50
    // for each degree of longitude between -66 and -125
    let range = CoordinateRange {
        min_longitude: -125,
        max_longitude: -66,
        min_latitude: 30,
        max_latitude: 48,
    };

    let a = angular_precessed_planets_in_range(radix_dt, target_dt, range, true, true, false);
    for key in a.as_ref().unwrap().keys() {
        println!("{}", key);
        let v = &a.as_ref().unwrap()[key];
        for angularity in v {
            let (degree, minute, _) = utils::decimal_to_dms(angularity.orb);
            println!(
                "{} {} {}\u{00B0}{}'",
                angularity.planet, angularity.framework, degree, minute
            );
        }
    }
    println!("Checked {} locations", a.unwrap().len());
    close_dll();
}
