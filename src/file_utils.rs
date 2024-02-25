// Nova Rust, a core library wrapping the Swiss Ephemeris library for astrological use.
// Copyright (C) 2024 Mike Verducci
/*This file is part of Nova Rust.

Nova Rust is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Nova Rust is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
*/

use std::fs::File;
use std::io::prelude::*;

pub fn create_file(title: &str) -> anyhow::Result<File> {
    let mut file = File::create(format!("{}.txt", title))?;
    Ok(file)
}

pub fn write_to_file(
    file: &mut File,
    text: String,
) -> anyhow::Result<()> {
    file.write_all(text.as_bytes())?;
    Ok(())
}