# The Nova (Rust) API
This is the beginning of a web API and core library wrapping the [Swiss Ephemeris](https://www.astro.com/swisseph/swephinfo_e.htm) astrological library.

For now, there are some command-line utilities of use to astrologers. Eventually, this will become an actual web backend for a variety of web apps.

###TODO list
- [ ] License as required by Swiss Ephemeris
- [ ] Windows installer
- [ ] MacOS installer
- [ ] Linux installers? (what resulting format to use; .rpm?)
- [ ] Use binary search instead of linear search to find solunar return datetimes
- [ ] Reduce cloning on vectors
- [ ] Clean up type system
- [ ] Better error messages instead of unwrapping