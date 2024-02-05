pub mod lib {

    use anyhow::{anyhow, Result};
    use chrono::{DateTime, Datelike, TimeZone, Timelike, Utc, NaiveDate};
    use chrono_tz::Tz;
    use libc::{c_char, c_int};
    use std::ffi::{c_double, CString};
    use std::ptr;

    const SIDEREALMODE: c_int = 64 * 1024;

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
    }
    
    pub struct CoordinateSystemValues {
        pub longitude : f64,
        pub latitude: f64,
        pub speed_in_longitude : f64,
        pub right_ascension : f64,
        pub prime_vertical_longitude : f64,
        pub azimuth: f64,
        pub meridian_longitude: f64, // This is experimental
    }
    
    pub struct Planet {
        pub name: String,
        pub body_number: i32,
        pub coordinates: CoordinateSystemValues,
    }
    
    #[derive(Debug)]
    pub struct Location {
        pub longitude: f64,
        pub latitude: f64,
        pub datetime: DateTime<Tz>,
        pub name: String,
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
    
    pub fn open_dll() {
        let path = "/home/mike/projects/nova_go_api/src/swe/ephemeris";
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

            if err_string.is_none() {
                Ok(float_array)
            } else {
                Err(anyhow!(err_string.unwrap()))
            }
        }
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
    
    fn error_string_from_buffer(err_buffer: &[c_char; 256]) -> Option<String> {
        for c in err_buffer {
            if *c != 0 {
                let err_cstring = unsafe { CString::from_raw(err_buffer.as_ptr() as *mut c_char) };
                let err_string = err_cstring.to_string_lossy().into_owned();
                return Some(err_string.trim_matches(char::from(0)).to_string());
            }
        }
        None
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

    fn calculate_lst(dt: DateTime<Utc>, decimal_longitude: f64) -> f64 {
        let julian_day_0_gmt = get_julian_day_0_gmt(dt);
    
        let universal_time =
            convert_dms_to_decimal(dt.hour() as f64, dt.minute() as f64, dt.second() as f64);
    
        let sidereal_time_at_midnight_julian_day = (julian_day_0_gmt - 2451545.0) / 36525.0;
    
        let greenwich_sidereal_time = 6.697374558
            + 2400.051336 * sidereal_time_at_midnight_julian_day
            + 0.000024862 * sidereal_time_at_midnight_julian_day.powi(2)
            + universal_time * 1.0027379093;
    
        let local_sidereal_time =
            (greenwich_sidereal_time + (decimal_longitude / 15.0)).rem_euclid(24.0);
    
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
    
    // fn try_date_with_timezone(dt: NaiveDateTime, tz_string: String) -> anyhow::Result<DateTime<Tz>> {
    //     let tz: Tz = tz_string.parse().ok().unwrap()?;
    //     let result = dt.and_local_timezone(tz);
    // }
    
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
                false => Ok(swe_array[0]),
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

}
