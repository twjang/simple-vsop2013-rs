#!/usr/bin/env python

import os

PATH_SELF = os.path.realpath(os.path.dirname(__file__))
PATH_SRC = os.path.join(PATH_SELF, '..', 'src')
PATH_BINS = os.path.join(PATH_SRC, 'bins')

def main():
    bodyvar2maxpower = {}

    for fname in os.listdir(PATH_BINS):
        if not fname.endswith('.bin'):
            continue
        _, bodyname, varid, t_power = fname.rsplit('.', 1)[0].split('_')
        key = (bodyname, varid)
        if not key in bodyvar2maxpower:
            bodyvar2maxpower[key] = int(t_power)
        else:
            bodyvar2maxpower[key] = max(int(t_power), bodyvar2maxpower[key])

    bodylst = [
        'mercury','venus','emb','mars','jupiter',
        'saturn','uranus','neptune','pluto'
    ]

    code = []
    code.append('// This file was generated using script/generate_code.py')
    code.append('')
    code.append('use super::utils::*;')
    code.append('use super::coeff;')
    code.append('')
    for bodyidx, body in enumerate(bodylst):
        code.append(f'/// Calculate elliptic variables for {body.upper()}. The frame is dynamical')
        code.append(f'/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore')
        code.append(f'/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.')
        code.append(f'///')
        code.append(f'/// # Return')
        code.append(f'///   * out.0: semi-major axis(au)')
        code.append(f'///   * out.1: mean logitude (rd)')
        code.append(f'///   * out.2: k = e * cos(pi) (rd)')
        code.append(f'///   * out.3: h = e * sin(pi) (rd)')
        code.append(f'///   * out.4: q = sin(i/2) * cos(omega) (rd)')
        code.append(f'///   * out.5: p = sin(i/2) * sin(omega) (rd)')
        code.append(f'///  where')
        code.append(f'///   *     e: eccentricity')
        code.append(f'///   *    pi: perihelion longitude')
        code.append(f'///   *     i: inclination')
        code.append(f'///   * omega: ascending node longitude')
        code.append(f'#[cfg(feature="{body}")]')
        code.append(f'pub fn {body}_elliptic_vars(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {{')
        code.append(f'    let tj = (t - coeff::CONST_T2000) / coeff::CONST_A1000;')
        code.append('')
        code.append(f'    let mut res: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);')
        code.append(f'    let mut tjp = 1.0;')
        code.append('')

        max_powers = [bodyvar2maxpower[(body, str(i))] for i in range(1, 7)]
        for cur_power in range(0, 21):
            cnt = 0

            cur_code_blk = []
            cur_code_blk.append(f'    // tj ^ {cur_power}')
            for varidx in range(0, 6):
                if max_powers[varidx] >= cur_power:
                    cur_code_blk.append(f'    res.{varidx} += accumulate(u64_as_f64_slice(&coeff::DATA_COEFF_{body.upper()}_{varidx+1}_{cur_power}), tj, tjp, tol);')
                    cnt += 1

            if cnt > 0:
                if cur_power > 0:
                    code.append('    tjp *= tj;')
                    code.append('')
                code.extend(cur_code_blk)

            if cnt == 0:
                break
        code.append('')
        code.append(f'    res.1 += coeff::CONST_FREQPLA[{bodyidx}] * tj;')
        code.append(f'    res.1 = mod_float(res.1, DPI);')
        code.append('')
        code.append('    res')
        code.append('}')
        code.append('')
        code.append(f'/// Calculate planetary rectangular coordinates for {body.upper()}. The frame is dynamical')
        code.append(f'/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore')
        code.append(f'/// terms whose `(cos coeff.^2 + sin coeff.^2) ^ 0.5` are smaller than `tol`.')
        code.append(f'///')
        code.append(f'/// # Return ')
        code.append(f'///  * out.0: X  (au)')
        code.append(f'///  * out.1: Y  (au)')
        code.append(f'///  * out.2: Z  (au)')
        code.append(f'///  * out.3: X\' (au/day)')
        code.append(f'///  * out.4: Y\' (au/day)')
        code.append(f'///  * out.5: Z\' (au/day)')
        code.append(f'#[cfg(feature="{body}")]')
        code.append(f'pub fn {body}_cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {{')
        code.append(f'   let ep = {body}_elliptic_vars(t, tol);')
        code.append(f'   elliptic_vars_to_cartesian(&ep, coeff::CONST_RGM[{bodyidx}])')
        code.append('}')
        code.append('')

    with open(os.path.join(PATH_SRC, 'generated.rs'), 'w') as f:
        f.write('\n'.join(code))
        f.write('\n')

if __name__ == '__main__':
    main()