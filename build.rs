use std::env;
use std::path::PathBuf;

fn main() {
    let lib_src_path = PathBuf::from("src").join("swe").join("dll");

    println!("cargo:rustc-link-search=native={}", lib_src_path.display());
    println!("cargo:rustc-link-lib=swe");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_src_path.display());
}
