use chrono::{NaiveDate, TimeZone};
use nova_rust::lib::{close_dll, create_sidereal_context, open_dll, Location};

fn main() {
    open_dll();
    let local = NaiveDate::from_ymd_opt(1989, 12, 20).unwrap().and_hms_opt(22, 20, 0).unwrap();
    let dt = chrono_tz::EST.from_local_datetime(&local).unwrap();

    let loc = Location {
        longitude:-74.1166,
        latitude: 40.97,
        datetime: dt,
        name: String::from("mv"),
    };
    let c = create_sidereal_context(dt, loc);
    close_dll();

    println!("{:?}", c);
}