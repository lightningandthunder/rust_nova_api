fn main() {
       // Specify the path where to find the native library
       println!("cargo:rustc-link-search=native=/home/mike/projects/nova-rust/src/swe/dll/");

       // Specify the native library name
       println!("cargo:rustc-link-lib=swe");

       println!("cargo:rustc-link-arg=-Wl,-rpath,/home/mike/projects/nova-rust/src/swe/dll/");
}