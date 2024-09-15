#!/usr/bin/env python

import os
import math
import cmath
import struct
import argparse
import urllib.request
import numpy as np

PATH_SELF = os.path.realpath(os.path.dirname(__file__))

np_dtype = np.float64

def mod_float(x, mod):
    return x - np.floor(x / mod) * mod

class VSOP2013:
    def __init__(self, ibody: int, fname: str):
        """
            Load the parameter files and prepare constants.

            @param {int} ibody specify parameters of which planet will be calculated.
!              1: Mercury
!              2: Venus
!              3: Earth-Moon barycenter
!              4: Mars
!              5: Jupiter
!              6: Saturn
!              7: Uranus
!              8: Neptune
!              9: Pluto
            @param {str} fname path to the parameter file.
        """
        self.ibody = ibody
        self.fname = fname
        self.ci0 =  [              # Mean Longitude J2000 (radian)
            0.4402608631669000e1,  # Mercury
            0.3176134461576000e1,  # Venus
            0.1753470369433000e1,  # Earth-Moon Barycenter
            0.6203500014141000e1,  # Mars 
            0.4091360003050000e1,  # Vesta
            0.1713740719173000e1,  # Iris
            0.5598641292287000e1,  # Bamberga
            0.2805136360408000e1,  # Ceres
            0.2326989734620000e1,  # Pallas
            0.5995461070350000e0,  # Jupiter
            0.8740185101070000e0,  # Saturn
            0.5481225395663000e1,  # Uranus
            0.5311897933164000e1,  # Neptune
            0.e0,                  # 
            5.19846640063e0,       # Moon (D)
            1.62790513602e0,       # Moon (F)
            2.35555563875e0        # Moon (l) 
        ]
        self.ci1= [                # Mean Motions in longitude (radian/cy)
            0.2608790314068555e5,  # Mercury
            0.1021328554743445e5,  # Venus
            0.6283075850353215e4,  # Earth-Moon Barycenter
            0.3340612434145457e4,  # Mars 
            0.1731170452721855e4,  # Vesta
            0.1704450855027201e4,  # Iris
            0.1428948917844273e4,  # Bamberga
            0.1364756513629990e4,  # Ceres
            0.1361923207632842e4,  # Pallas
            0.5296909615623250e3,  # Jupiter
            0.2132990861084880e3,  # Saturn
            0.7478165903077800e2,  # Uranus
            0.3813297222612500e2,  # Neptune
            0.3595362285049309e0,  # Pluto (Mu from TOP2013)
              77713.7714481804e0,  # Moon (D)
              84334.6615717837e0,  # Moon (F)
              83286.9142477147e0   # Moon (l)
        ]

        self.freqpla=[ 	           # Planetary frequency in longitude
            0.2608790314068555e5,  # Mecrcury
            0.1021328554743445e5,  # Venus
            0.6283075850353215e4,  # Earth-Moon Barycenter
            0.3340612434145457e4,  # Mars 
            0.5296909615623250e3,  # Jupiter
            0.2132990861084880e3,  # Saturn
            0.7478165903077800e2,  # Uranus
            0.3813297222612500e2,  # Neptune
            0.2533566020437000e2,  # Pluto
        ]
#
# --- Rotation Matrix : Ecliptic -> Equator ----------------------------
#
        dgrad = math.pi / 180.0
        sdrad = math.pi / 180.0 / 3600.0

        # eps = 23 26.0' 21.41136''
        # phi = -0.05188''
        eps  = (23.e0 + 26.e0/60.e0 + 21.41136e0/3600.e0)*dgrad
        phi  = -0.05188e0 * sdrad

        ceps = math.cos(eps)
        seps = math.sin(eps)
        cphi = math.cos(phi)
        sphi = math.sin(phi)

        rot = np.array([
            cphi, -sphi*ceps, sphi*seps,
            sphi, cphi*ceps, -cphi*seps,
            0.0 , seps, ceps
        ], dtype=np_dtype).reshape(3, 3)
        self.rot = rot

# --- Masses system (INPOP10A) -----------------------------------------
        self.gmp = [
            4.9125474514508118699e-11,
            7.2434524861627027000e-10,
            8.9970116036316091182e-10,
            9.5495351057792580598e-11,
            2.8253458420837780000e-07,
            8.4597151856806587398e-08,
            1.2920249167819693900e-08,
            1.5243589007842762800e-08,
            2.1886997654259696800e-12
        ]
        self.gmsol = 2.9591220836841438269e-04
# ----------------------------------------------------------------------


        self.a1000 = 365250.e0
        self.t2000 = 2451545.e0

        self.maxterm=351000
        self.limit = np.zeros([7, 21], dtype=np.int32) 
        self.iphi = np.zeros([17, self.maxterm], dtype=np_dtype)
        self.cc = np.zeros([self.maxterm], dtype=np_dtype)
        self.ss = np.zeros([self.maxterm], dtype=np_dtype)
        self.aa = np.zeros([self.maxterm], dtype=np_dtype)
        self.bb = np.zeros([self.maxterm], dtype=np_dtype)
        
        ip = 0
        iv = 0
        it = 0
        nt = 0

        nn = 0
        remaining_line = 0
        f = open(fname, 'r')

        while True:
            line = f.readline()
            if line == '': return

            ip = int(line[9:12].strip())
            iv = int(line[12:15].strip())
            it = int(line[15:18].strip())
            nt = int(line[18:25].strip())

            assert ibody == ip, f'file planet number({ip}) is not equal to fn arg({ibody})'
            self.limit[iv, it] = nt
            
            for n in range(1, nt+1):
                line = f.readline()
                if line == '': return
                if iv == 2 and it == 1 and n == 1:
                    self.limit[iv, it] = nt-1
                    continue

                values = []
                values.extend(int(line[i:i+3]) for i in range(6, 18, 3))  # 4 3-digit integers
                values.extend(int(line[i:i+3]) for i in range(19, 34, 3))  # 5 3-digit integers
                values.extend(int(line[i:i+4]) for i in range(35, 51, 4))  # 4 4-digit integers
                values.append(int(line[52:58]))  # 1 6-digit integer
                values.extend(int(line[i:i+3]) for i in range(59, 68, 3))  # 3 3-digit integers
                values.append(line[68:92])
                values.append(line[92:116])  # 2 24-character strings
                values[17]=float(values[17][:20]+'e'+values[17][21:])
                values[18]=float(values[18][:20]+'e'+values[18][21:])

                for idx in range(17):
                    self.iphi[idx, nn] = values[idx]
                    self.aa[nn] += self.iphi[idx, nn] * self.ci0[idx]
                    self.bb[nn] += self.iphi[idx, nn] * self.ci1[idx]

                self.ss[nn] = values[17]
                self.cc[nn] = values[18]
                
                nn += 1


    def calculate(self, tdj: float)->np.ndarray:
        """
            Calculate elliptic variables as described below. The frame is dynamical equinox and ecliptic J2000.
            The time scale is TDB(Temps Dynamique Barycentrique)

            @param {float} tdj julian date in dynamical time TDB from J2000
            @returns {np.ndarray} elliptic variable of the given planet.
                [0]: semi-major axis (au)
                [1]: mean longitude (rd)
                [2]: k = e*cos(pi) (rd)
                [3]: h = e*sin(pi) (rd)
                [4]: q = sin(i/2)*cos(omega) (rd)
                [5]: p = sin(i/2)*sin(omega) (rd)
                e:     eccentricity
                pi:    perihelion longitude
                i:     inclination
                omega: ascending node longitude
        """
        tj = (tdj - self.t2000) / self.a1000
        t = np.power(tj, np.arange(0, 21).astype(np_dtype)).astype(np_dtype)
        t[0] = 1.0

        res = np.zeros([6], dtype=np_dtype)
        nn = 0
        dpy = np.pi * 2.0

        for idx_output in range(1, 7):
            cnt = 0
            for idx_t_power in range(21):
                if self.limit[idx_output, idx_t_power] == 0:
                    continue
                for n in range(self.limit[idx_output, idx_t_power].item()):
                    aa=self.aa[nn]
                    bb=self.bb[nn]
                    arg = aa + bb * t[1]
                    arg = mod_float(arg, dpy)
                    res[idx_output-1] += t[idx_t_power] * ( 
                                          self.ss[nn] * np.sin(arg) + 
                                          self.cc[nn] * np.cos(arg)
                                          )
                    nn += 1
        
        xl = res[1] + self.freqpla[self.ibody-1] * tj
        xl = mod_float(xl, np.pi * 2.0)
        res[1] = xl

        return res

    def ellxyz(self, v: np.ndarray)->np.ndarray:
        """
            Compute planetary rectangular coordinates from elliptic variables.

            @param {np.ndarray} v elliptic variables returned from calculate() function
            @returns {np.ndarray} rectangular coordinates
!              w[i],i=0..3 Positions  X, Y, Z (au)
!              w[i],i=3..6 Velocities X',Y',Z'(au/d)
        """
#
# --- Initialization ---------------------------------------------------
#
        rgm = (self.gmp[self.ibody-1] + self.gmsol) ** 0.5
        xa, xl, xk, xh, xq, xp = v.reshape(-1).tolist()
#
# --- Computation ------------------------------------------------------
#
        xfi = (1.e0 - xk*xk - xh*xh) ** 0.5
        xki = (1.e0 - xq*xq - xp*xp) ** 0.5
        u  =1.e0/(1.e0+xfi)

        z = complex(xk, xh)
        ex = abs(z)
        ex2 = ex*ex
        ex3 = ex2*ex
        z1 = z.conjugate()

        gl = mod_float(xl, 2 * math.pi)
        gm = gl - math.atan2(xh, xk)
        e = gl + (ex-0.125*ex3)*math.sin(gm) + 0.5*ex2*math.sin(2.0*gm) + 0.375*ex3*math.sin(3.0*gm)
        
        while True:
            z2 = complex(0.0, e)
            zteta = cmath.exp(z2)
            z3 = z1 * zteta
            dl = gl - e + z3.imag
            rsa = 1.e0 - z3.real
            e = e + dl/rsa
            if abs(dl) < 1e-15: break
        
        z1 = u * z * z3.imag 
        z2 = complex(z1.imag, -z1.real)
        zto = (-z + zteta + z2) / rsa
        xcw = zto.real
        xsw = zto.imag
        xm = xp * xcw - xq * xsw
        xr = xa * rsa

        w0 = xr * (xcw - 2.0 * xp * xm)
        w1 = xr * (xsw + 2.0 * xq * xm)
        w2 = -2.0 * xr * xki * xm

        xms = xa * (xh+xsw) / xfi
        xmc = xa * (xk+xcw) / xfi
        xn = rgm / (xa ** 1.5)

        w3 = xn * ((2.0*xp*xp-1.0)*xms + 2.0*xp*xq*xmc)
        w4 = xn * ((1.0-2.0*xq*xq)*xmc - 2.0*xp*xq*xms)
        w5 = 2.0*xn*xki*(xp*xms + xq*xmc)

        res = np.array([w0, w1, w2, w3, w4, w5], dtype=np_dtype)
        return res
   
    
    def elliptic_to_equator(self, res: np.ndarray)->np.ndarray:
        """
            Convert ellxyz output(elliptic coordinates) to equator coordinates

            @param {np.ndarray} res outputs of ellxyz()
            @returns {np.ndarray} 
!              w[i],i=0..3 Positions  X, Y, Z (au)
!              w[i],i=3..6 Velocities X',Y',Z'(au/d)
        """
        new_res1 = np.matmul(self.rot, res[0:3])
        new_res2 = np.matmul(self.rot, res[3:6])
        res = np.concatenate([new_res1, new_res2])
        return res
    
    def export(self)->dict:
        nn = 0
        coeff = {}

        for idx_output in range(1, 7):
            for idx_t_power in range(21):
                cur_coeff = []
                for n in range(self.limit[idx_output, idx_t_power].item()):
                    aa = self.aa[nn]
                    bb = self.bb[nn]
                    ss = self.ss[nn]
                    cc = self.cc[nn]
                    cur_coeff.append((ss, cc, aa, bb))
                    nn += 1
                cur_coeff.sort(key=lambda e: (e[0] * e[0] + e[1] * e[1])**0.5)
                coeff[(idx_output, idx_t_power)] = cur_coeff


        return {
            'coeff': coeff,
            'common': {
                't2000': self.t2000,
                'a1000': self.a1000,
                'freqpla': self.freqpla,
                'rgm': [(gmp + self.gmsol)**0.5 for gmp in self.gmp],
                'rot': self.rot.reshape(-1).tolist(),
            }
        }


def main_orig(args: argparse.Namespace):
    f = open('VSOP2013.py.out', 'w')
    body = [
        'MERCURY','VENUS','EMB','MARS','JUPITER',
        'SATURN','URANUS','NEPTUNE','PLUTO'
    ]

    fname = [
        'VSOP2013p1.dat','VSOP2013p2.dat',
        'VSOP2013p3.dat','VSOP2013p4.dat',
        'VSOP2013p5.dat','VSOP2013p6.dat',
        'VSOP2013p7.dat','VSOP2013p8.dat',
        'VSOP2013p9.dat'
    ]

    step = 4000.0
    ndat = 11
    t0   = 2411545.0e0   ###  1890 June 26 12h

    for ip in range(1, 10):
        print(f'  *** {body[ip-1]}')
        vsop2013 = VSOP2013(ip, os.path.join(args.data, fname[ip-1]))
        f.write(f"""
  PLANETARY EPHEMERIS VSOP2013  {body[ip-1]}

  1/ Elliptic   Elements:    a (au), lambda (radian), k, h, q, p        - Dynamical Frame J2000
  2/ Ecliptic   Heliocentric Coordinates:  X,Y,Z (au)  X',Y',Z' (au/d)  - Dynamical Frame J2000
  3/ Equatorial Heliocentric Coordinates:  X,Y,Z (au)  X',Y',Z' (au/d)  - ICRS Frame J2000

""")
        for n in range(1, ndat+1):
            t = t0 + (n-1) * step
            el = vsop2013.calculate(t)
            r1 = vsop2013.ellxyz(el)
            r2 = vsop2013.elliptic_to_equator(r1)
            f.write(f'  Julian Date JD{t:10.1f} (TDB)\n')

            e = [f'{e:16.10f}' for e in el.tolist()]
            f.write(''.join(e) + '\n')

            e = [f'{e:16.10f}' for e in r1.tolist()]
            f.write(''.join(e) + '\n')

            e = [f'{e:16.10f}' for e in r2.tolist()]
            f.write(''.join(e) + '\n')
            f.flush()



def main_export(args: argparse.Namespace):
    path_export = args.output
    if path_export is None:
        raise ValueError('--output is required')
    
    body = [
        'MERCURY','VENUS','EMB','MARS','JUPITER',
        'SATURN','URANUS','NEPTUNE','PLUTO'
    ]

    fname = [
        'VSOP2013p1.dat','VSOP2013p2.dat',
        'VSOP2013p3.dat','VSOP2013p4.dat',
        'VSOP2013p5.dat','VSOP2013p6.dat',
        'VSOP2013p7.dat','VSOP2013p8.dat',
        'VSOP2013p9.dat'
    ]

    os.makedirs(path_export, exist_ok=True)
    os.makedirs(os.path.join(path_export, 'bins'), exist_ok=True)

    export = {
        'coeff': dict(),
    }

    for ip in range(1, 10):
        vsop2013 = VSOP2013(ip, os.path.join(args.data, fname[ip-1]))
        cur_export = vsop2013.export()
        export.update(cur_export['common'])
        export['coeff'][body[ip-1]] = cur_export['coeff']

    max_t_power = {}
    for ip in range(1, 10):
        cur_body = body[ip-1]
        for (idx_output, idx_t_power), v in export['coeff'][cur_body].items():
            if not (ip, idx_output) in max_t_power:
                max_t_power[(ip, idx_output)] = 0
            if len(v) > 0:
                max_t_power[(ip, idx_output)] = max(idx_t_power, max_t_power[(ip, idx_output)])

    for ip in range(1, 10):
        cur_body = body[ip-1]
        for (idx_output, idx_t_power), v in export['coeff'][cur_body].items():
            if idx_t_power > max_t_power[(ip, idx_output)]: continue
            data = []
            num_data = len(v)

            for idx_start in range(0, num_data, 4):
                idx_end = min(idx_start + 4, num_data)
                cur_batch = v[idx_start: idx_end].copy()
                if len(cur_batch) < 4:
                    cur_batch.extend([(0.0, 0.0, 0.0, 0.0)] * (4 - len(cur_batch)))

                # we will transform 
                # cur_batch = [
                #  (s0, c0, a0, b0), (s1, c1, a1, b1), (s2, c2, a2, b2), (s3, c3, a3, b3),   
                # ] to 
                # cur_batch = [
                #  (s0, s1, s2, s3), (c0, c1, c2, c3), (a0, a1, a2, a3), (b0, b1, b2, b3),   
                # ]
                # to enable simd 
                new_cur_batch = []
                for i in range(4):
                    for j in range(4):
                        new_cur_batch.append(cur_batch[j][i])
                sz = (cur_batch[0][0]*cur_batch[0][0] + cur_batch[0][1]*cur_batch[0][1]) ** 0.5
                data.append(sz)
                data.extend(new_cur_batch)
                

            fname = f'coeff_{cur_body.lower()}_{idx_output}_{idx_t_power:02d}.bin'
            with open(os.path.join(path_export, 'bins', fname), 'wb') as f:
                for e in data:
                    f.write(struct.pack('d', e))

    code = []
    code.append('// This file was generated using script/vsop2013.py')
    code.append('')
    code.append('use include_bytes_plus::include_bytes;')
    code.append('')
    code.append('// All elements of coeff. arrays are f64 and follows following format:')
    code.append('// [sz=|(s0, c0)|, (s0, s1, s2, s3), (c0, c1, c2, c3), (a0, a1, a2, a3), (b0, b1, b2, b3)] ...repeated... ') 
    code.append('')
    for ip in range(1, 10):
        cur_body = body[ip-1]
        for (idx_output, idx_t_power), v in export['coeff'][cur_body].items():
            if idx_t_power > max_t_power[(ip, idx_output)]: continue
            n_data = ((len(v) + 3) // 4) * 17

            fname = f'coeff_{cur_body.lower()}_{idx_output}_{idx_t_power:02d}.bin'
            code.append(f'#[cfg(feature="{cur_body.lower()}")]')
            code.append(f'pub const DATA_COEFF_{cur_body.upper()}_{idx_output}_{idx_t_power}: [u64; {n_data}] = include_bytes!("src/bins/{fname}" as u64le);')
    
    for k, v in export.items():
        if k == 'coeff':
            continue
        if isinstance(v, float):
            code.append(f'pub const CONST_{k.upper()}: f64 = {v};');
        elif isinstance(v, list):
            code.append(f'pub const CONST_{k.upper()}: [f64;{len(v)}] = [');
            v = [f'{e:27.20e}' for e in v]
            for idx_start in range(0, len(v), 3):
                idx_end = min(idx_start + 3, len(v))
                line = ' ' * 2
                line += ', '.join(v[idx_start:idx_end])
                if idx_end < len(v):
                    line += ','
                code.append(line)
            code.append('];')
    
    with open(os.path.join(path_export, 'coeff.rs'), 'w') as f:
        f.write('\n'.join(code)) 

def main_download(args: argparse.Namespace):
    base_url = 'https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/'
    fnames = [
        'README.pdf',
        'VSOP2013-secular.dat',
        'VSOP2013.ctl',
        'VSOP2013.f',
        'VSOP2013p1.dat',
        'VSOP2013p2.dat',
        'VSOP2013p3.dat',
        'VSOP2013p4.dat',
        'VSOP2013p5.dat',
        'VSOP2013p6.dat',
        'VSOP2013p7.dat',
        'VSOP2013p8.dat',
        'VSOP2013p9.dat',
    ]
    for fname in fnames:
        url = base_url + fname
        path_dst = os.path.join(args.data, fname)
        os.makedirs(os.path.dirname(path_dst), exist_ok=True)

        if os.path.exists(path_dst):
            print(f'! skipping {url}: {path_dst} exists')
        else:
            print(f'* downloading {url} -> {path_dst}')
            urllib.request.urlretrieve(url, path_dst)
    print(f'* done')



def main_check(args: argparse.Namespace):
    import subprocess

    main_download(args)
    subprocess.call([
        'gfortran', 'VSOP2013.f', '-o', 'VSOP2013', '-O3'
    ], cwd=args.data)
    subprocess.call([
        './VSOP2013'
    ], cwd=args.data)

    main_orig(args)
    
    with open('VSOP2013.py.out') as f:
        res_ported = [x.strip() for x in f.readlines() if x.strip() != '']
    with open(os.path.join(args.data, 'VSOP2013.out')) as f:
        res_orig = [x.strip() for x in f.readlines() if x.strip() != '']
    
    assert len(res_ported) == len(res_orig)
    for line_ported, line_orig in zip(res_ported, res_orig):
        assert line_ported == line_orig

    print('* check passed!')



 
def parse_args()->argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['main', 'export', 'download', 'check'])
    parser.add_argument('-d,--data', help='path to the directory containing VSOP2013p{1..9} files', 
                        dest='data', 
                        default=os.path.join(PATH_SELF, 'data'))
    parser.add_argument('-o,--output', help='path to output', dest='output')
    args = parser.parse_args()
    return args

           
if __name__ == '__main__':
    args = parse_args()
    mode = args.mode

    if mode == 'main':
        main_orig(args)
    elif mode == 'export':
        main_export(args)
    elif mode == 'download':
        main_download(args)
    elif mode == 'check':
        main_check(args)
    else:
        raise ValueError('unsupported mode {args.mode}')