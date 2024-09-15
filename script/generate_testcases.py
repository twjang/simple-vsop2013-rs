#!/usr/bin/env python
"""
Running `vsop2013.py main` generates `VSOP2013.py.out`. This script transforms the output file 
to test cases.
"""
import os

PATH_SRC = os.path.realpath(os.path.dirname(__file__))

def main():
    with open(os.path.join(PATH_SRC, 'VSOP2013.py.out')) as f:
        lines = f.readlines()

    cases: dict[str, dict[float, list]] = {}
    cur_planet = None

    for line in lines:
        if line.find('PLANETARY EPHEMERIS') > -1:
            cur_planet = line.strip().split()[-1]
            if cur_planet not in cases:
                cases[cur_planet] = {}
            continue
        elif line.find('(TDB)') > -1:
            cur_t = float(line.strip().split()[-2])
            if cur_t not in cases[cur_planet]:
                cases[cur_planet][cur_t] = []
            continue
        elif len(line.strip()) > 0:
            is_fp = False
            try:
                values = [float(e) for e in line.strip().split()]
                is_fp = len(values) == 6
            except ValueError:
                values = None
            if is_fp:
                cases[cur_planet][cur_t].append(line.strip().split())

    case_code = ''

    for planet in cases.keys():
        dates = sorted(cases[planet].keys())
        case_code +=f'    #[cfg(feature="{planet.lower()}")]\n'
        case_code +=f'    #[test]\n'
        case_code +=f'    fn check_{planet.lower()}() {{\n'
        case_code += '        let cases = vec![\n'
        for cur_date in dates:
            case_code += f'            ({cur_date:.1f}, \n'
            case_code += '                (' + ', '.join(cases[planet][cur_date][0]) + '),\n'
            case_code += '                (' + ', '.join(cases[planet][cur_date][1]) + '),\n'
            case_code += '            ),\n'
        case_code += '        ];\n'
        case_code += f"""
        let mut failcnt = 0;
        for (t, ep_expected, xyz_expected) in cases.iter() {{
            let ep_actual = {planet.lower()}_elliptic_vars(*t, 0.0);
            let xyz_actual = {planet.lower()}_cartesian(*t, 0.0);

            let mut failed = false;
            if !is_matched(ep_expected, &ep_actual) {{
                failed = true;
            }}
            if !is_matched(xyz_expected, &xyz_actual) {{
                failed = true;
            }}
            if !failed {{
                println!("{planet}(t={{t}}) PASSED")
            }} else {{
                failcnt += 1; 
                println!("{planet}(t={{t}}) FAILED")
            }}
        }}

        assert_eq!(failcnt, 0);\n"""
        case_code += '    }\n'

    print(case_code)

if __name__ == '__main__':
    main()