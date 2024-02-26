use std::env;
use std::path::PathBuf;

fn main() {
    let target = env::var("TARGET").unwrap();
    let current_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    let lib_dir = PathBuf::from(&current_dir).join("src/swe/dll/");
    let lib_dir_str = lib_dir.to_str().unwrap();

    println!("cargo:rustc-link-search=native={}", lib_dir_str);

    println!("cargo:rustc-link-lib=swe");

    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_dir_str);
}