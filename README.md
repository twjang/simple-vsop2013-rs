# simple-vsop2013 

`VSOP2013` provides accurate positions for the eight major planets of the Solar System (Mercury to Neptune) over a time span of several thousand years. This Rust library is a port of the original VSOP2013 Fortran implementation, which can be found at https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/. 

Several optimizations have been applied:
* SIMD 
* Ignoring small terms using the `tol` parameter
* Including only necessary coefficients

Like the original implementation, this library operates based on the frame defined by the dynamical equinox and ecliptic J2000. The time scale of `t` is TDB (Barycentric Dynamical Time).

## Usage

### Basic Usage

Here's a basic example of how to use the library:

```rust
use simple_vsop2013;

fn main() {
    const J2000: f64 = 2451545.0; // 2000-01-01, 12:00:00 TT
    let coords = simple_vsop2013::emb_cartesian(J2000, 0.0);
    
    println!("Earth's position on 2000-01-01:");
    println!("X : {} AU", coords.0);
    println!("Y : {} AU", coords.1);
    println!("Z : {} AU", coords.2);
    println!("X': {} AU/d", coords.3);
    println!("Y': {} AU/d", coords.4);
    println!("Z': {} AU/d", coords.5);
}
```

### Features

You can selectively include coefficients for each planet through **features**. 
For example, if you enable the 'venus' feature, only Venus-related functions and Venus-related coefficients are included during library compilation. If you enable both 'venus' and 'emb', functions related to both planets are included during compilation. 

The following features are supported:
- mercury
- venus
- emb
- mars
- jupiter
- saturn
- uranus
- neptune
- pluto

### Scripts

* `generate_code.py`: Generates `{planet}_cartesian()` and `{planet}_elliptic_vars` functions from the `.bin` files in the `src/bins` directory.
* `vsop2013.py`: Supports the following features:
  * Downloads the official VSOP2013 implementation and coefficients from FTP.
  * Replicates the functionality of `VSOP2013.f`
  * Generates `.bin` files. You may need to modify this file if you want to reduce the bin file size.
* `generate_testcases.py`: Reads `VSOP2013.py.out` (the output of `vsop2013.py`) and prints test functions.


## Acknowledgements
This implementation is based on the VSOP2013 theory developed by G. Francou and J.-L. Simon at the Institut de Mécanique Céleste et de Calcul des Éphémérides (IMCCE).