// this file was generated using script/generate_code.py

use super::utils::*;
use super::coeff;

/// Calculate elliptic variables for MERCURY. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="mercury")]
pub fn mercury_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_8), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_5_8), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_6_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_1_10), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_3_10), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MERCURY_4_10), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[0] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for MERCURY. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="mercury")]
pub fn mercury_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = mercury_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[0])
}

/// Calculate elliptic variables for VENUS. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="venus")]
pub fn venus_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_8), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_5_8), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_6_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_4_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_VENUS_3_10), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[1] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for VENUS. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="venus")]
pub fn venus_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = venus_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[1])
}

/// Calculate elliptic variables for EMB. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="emb")]
pub fn emb_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_4_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_EMB_2_10), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[2] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for EMB. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="emb")]
pub fn emb_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = emb_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[2])
}

/// Calculate elliptic variables for MARS. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="mars")]
pub fn mars_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_5_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_4_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_MARS_2_9), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[3] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for MARS. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="mars")]
pub fn mars_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = mars_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[3])
}

/// Calculate elliptic variables for JUPITER. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="jupiter")]
pub fn jupiter_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_8), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_5_8), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_6_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_10), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_3_10), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_4_10), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 11
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_11), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_11), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 12
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_12), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_12), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 13
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_13), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_13), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 14
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_14), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_14), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 15
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_15), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_15), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 16
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_16), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_16), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 17
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_17), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_17), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 18
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_18), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_2_18), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 19
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_19), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 20
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_JUPITER_1_20), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[4] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for JUPITER. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="jupiter")]
pub fn jupiter_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = jupiter_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[4])
}

/// Calculate elliptic variables for SATURN. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="saturn")]
pub fn saturn_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_8), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_5_8), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_6_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_10), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_3_10), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_4_10), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 11
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_11), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_11), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 12
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_12), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_12), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 13
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_13), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_13), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 14
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_14), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_14), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 15
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_15), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_15), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 16
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_16), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_16), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 17
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_17), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_17), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 18
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_18), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_2_18), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 19
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_19), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 20
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_SATURN_1_20), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[5] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for SATURN. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="saturn")]
pub fn saturn_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = saturn_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[5])
}

/// Calculate elliptic variables for URANUS. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="uranus")]
pub fn uranus_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_1_10), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_3_10), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_URANUS_4_10), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[6] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for URANUS. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="uranus")]
pub fn uranus_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = uranus_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[6])
}

/// Calculate elliptic variables for NEPTUNE. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="neptune")]
pub fn neptune_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_8), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_5_8), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_9), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_10), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_3_10), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_4_10), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_6_10), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 11
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_NEPTUNE_1_11), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[7] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for NEPTUNE. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="neptune")]
pub fn neptune_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = neptune_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[7])
}

/// Calculate elliptic variables for PLUTO. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return
///   * out.0: semi-major axis(au)
///   * out.1: mean logitude (rd)
///   * out.2: k = e * cos(pi) (rd)
///   * out.3: h = e * sin(pi) (rd)
///   * out.4: q = sin(i/2) * cos(omega) (rd)
///   * out.5: p = sin(i/2) * sin(omega) (rd)
///  where
///   *     e: eccentricity
///   *    pi: perihelion longitude
///   *     i: inclination
///   * omega: ascending node longitude
#[cfg(feature="pluto")]
pub fn pluto_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;

    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut tjp = 1.0;

    // tj ^ 0
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_0), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_0), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_0), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_0), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_0), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_0), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 1
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_1), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_1), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_1), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_1), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_1), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_1), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 2
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_2), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_2), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_2), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_2), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_2), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_2), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 3
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_3), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_3), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_3), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_3), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_3), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_3), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 4
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_4), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_4), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_4), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_4), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_4), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_4), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 5
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_5), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_5), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_5), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_5), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_5), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_5), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 6
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_6), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_6), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_6), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_6), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_6), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_6), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 7
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_7), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_7), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_7), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_7), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_7), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_7), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 8
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_8), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_8), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_8), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_8), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_8), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_8), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 9
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_9), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_9), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_9), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_9), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_9), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_9), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 10
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_10), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_10), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_10), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_10), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_10), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_10), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 11
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_11), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_11), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_11), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_11), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_11), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_11), tj, tjp, tol);
    tjp *= tj;

    // tj ^ 12
    res.0 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_1_12), tj, tjp, tol);
    res.1 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_2_12), tj, tjp, tol);
    res.2 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_3_12), tj, tjp, tol);
    res.3 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_4_12), tj, tjp, tol);
    res.4 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_5_12), tj, tjp, tol);
    res.5 += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_PLUTO_6_12), tj, tjp, tol);

    res.1 += coeff::CONST_FREQPLA[8] * tj;
    res.1 = mod_float(res.1, DPI);

    res
}

/// Calculate planetary rectangular coordinates for PLUTO. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (au)
///  * out.1: Y  (au)
///  * out.2: Z  (au)
///  * out.3: X' (au/day)
///  * out.4: Y' (au/day)
///  * out.5: Z' (au/day)
#[cfg(feature="pluto")]
pub fn pluto_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
   let ep = pluto_elliptic_vars(t, tol);
   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[8])
}

