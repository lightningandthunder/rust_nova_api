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