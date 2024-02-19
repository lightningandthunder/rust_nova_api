use chrono::{DateTime, NaiveDate, NaiveDateTime, TimeZone, Utc};
use chrono_tz::Tz;
use console::Style;
use dialoguer::{
    theme::{ColorfulTheme, Theme},
    Confirm, Input, MultiSelect, Select,
};
use dotenv::dotenv;
use nova_rust::{
    angular_precessed_planets_in_range, close_dll, create_sidereal_context,
    find_next_solunar_utc_dt, open_dll, PLANET_NAMES, PLANET_TO_INT
};

use nova_rust::{spacetime_utils::{geocode, local_to_utc, naive_to_local_tz, string_to_naive_date, string_to_naive_datetime}, structs::CoordinateRange};

#[derive(Debug)]
struct InputConfig {
    title: String,
    local_birth_dt: DateTime<Tz>,
    birth_coordinates: (f64, f64),
    target_planet_name: String,
    harmonic: i32,
    // Will need to add 0 H/M/S to this
    solunar_dt: DateTime<Tz>,
    allowed_body_numbers: Vec<i32>,
    // Will need to convert to f64 if no decimal point is present
    orb: f64,
}

fn get_hand_entered_coordinates(theme: &ColorfulTheme) -> anyhow::Result<(f64, f64)> {
    let coordinate_input: String = Input::with_theme(theme)
    .with_prompt("Please enter longitude and latitude as signed decimals separated by a space; for example, -145.0091 34.0182")
    .interact()?;

    let raw: Vec<f64> = coordinate_input
        .trim()
        .split_whitespace()
        .into_iter()
        .map(|s| s.parse::<f64>().unwrap())
        .collect();

    Ok((raw[0], raw[1]))
}

fn get_natal_precession_input() -> anyhow::Result<Option<InputConfig>, Box<dyn std::error::Error>> {
    // Ask for title
    let theme = ColorfulTheme {
        values_style: Style::new().blue().dim(),
        ..ColorfulTheme::default()
    };
    println!("Welcome to the Nova Natal Planet Finder tool!");

    let title: String = Input::with_theme(&theme)
        .with_prompt("Report title")
        .interact()?;

    // Ask for birth date
    let birth_datetime_raw: String = Input::with_theme(&theme)
        .with_prompt("Please enter your birth date and time in this format: MM/DD/YYYY HH:MM am")
        .interact()?;

    let birth_time = string_to_naive_datetime(birth_datetime_raw)?;

    // Ask for birth location
    let birth_location: String = Input::with_theme(&theme)
        .with_prompt("Please enter the name of your birth location; for example: Jackson, MS, or: San Francisco, California.")
        .interact()?;

    let geocoded_coords = geocode(&birth_location)?;
    let birth_coordinates = match geocoded_coords {
        Some(coords) => {
            let prompt = format!("These coordinates were found for your birth location: {:?}. Enter Y to accept or N to enter your own coordinates.", coords);
            let use_geocoded_coordinates =
                Confirm::with_theme(&theme).with_prompt(prompt).interact()?;

            match use_geocoded_coordinates {
                true => coords,
                false => get_hand_entered_coordinates(&theme)?,
            }
        }
        None => get_hand_entered_coordinates(&theme)?,
    };

    let local_birth_dt = naive_to_local_tz(birth_time, birth_coordinates.0, birth_coordinates.1)?;

    // Ask for target planet
    let target_planet_int: usize = Select::with_theme(&theme)
        .with_prompt("Please select a planet to find a return for")
        .default(0)
        .item("Sun")
        .item("Moon")
        .interact()?;

    let target_planet_name = match target_planet_int {
        0 => "Sun".to_string(),
        _ => "Moon".to_string(),
    };

    // Ask for harmonic
    let harmonic_raw: String = Input::with_theme(&theme)
    .with_prompt("Please enter a harmonic to search for; for example, 1 for a full Solar or Lunar Return, 2 for a Demi, etc")
    .interact()?;

    let harmonic: i32 = harmonic_raw.parse().unwrap();

    // Ask for date after which to find the next return
    let target_date_raw: String = Input::with_theme(&theme)
        .with_prompt("Please enter the target date to search after in this format: MM/DD/YYYY")
        .interact()?;

    let solunar_start_date = string_to_naive_date(target_date_raw)?;

    let solunar_dt = naive_to_local_tz(solunar_start_date.and_hms_opt(0, 0, 0).unwrap(), birth_coordinates.0, birth_coordinates.1)?;

    // Select allowed planets
    let planet_defaults = [
        true,  // Sun
        false, // Moon
        false, // Mercury
        true,  // Venus
        false, // Mars
        true,  // Jupiter
        false, // Saturn
        true,  // Uranus
        false, // Neptune
        false, // Pluto,
        false, // Quaoar
        false, // Sedna
        false, // Orcus
        false, // Haumea
        false, // Eris
        false, // Makemake
        false, // Gonggong
    ];

    let selections = MultiSelect::with_theme(&theme)
        .with_prompt("Please select natal planets to include in your search. (Use up and down arrows to select planets, and spacebar to toggle them. Press enter when finished.)")
        .items(&PLANET_NAMES)
        .defaults(&planet_defaults)
        .interact()?;

    let allowed_body_numbers: Vec<i32> = selections
        .into_iter()
        .map(|planet_name_index| {
            let planet_name = PLANET_NAMES[planet_name_index];

            PLANET_TO_INT[planet_name]
        })
        .collect();

    // Select orb
    let orb_raw: String = Input::with_theme(&theme)
        .with_prompt(
            "Please enter the orb to allow for angular planets as a decimal, e.g. 1.0 or 2.5",
        )
        .default("3.0".to_string())
        .interact()?;

    let orb: f64 = orb_raw.parse().unwrap();

    println!("{:?}", String::from("3").parse::<f64>());

    Ok(Some(InputConfig {
        title,
        local_birth_dt,
        birth_coordinates,
        target_planet_name,
        harmonic,
        solunar_dt,
        allowed_body_numbers,
        orb,
    }))
}

fn main() {
    dotenv().ok();

    open_dll();

    let range = CoordinateRange {
        min_longitude: -125,
        max_longitude: -66,
        min_latitude: 30,
        max_latitude: 48,
    };

    let config = get_natal_precession_input().unwrap().unwrap();
    println!("{:?}", config);

    let r = angular_precessed_planets_in_range(config.local_birth_dt, config.solunar_dt, range, config.allowed_body_numbers).unwrap();
    println!("{:?}", r);

    close_dll();
}
