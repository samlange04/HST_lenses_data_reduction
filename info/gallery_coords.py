import re

# Source: Shu et al. 2016 (doi:10.3847/1538-4357/833/2/264), Table 1, column 1
# ("Selected Properties of the BELLS GALLERY Sample")
# Raw format: HHMMSS.ss[+/-]DDMMSS.s (J2000)
# Only E-S-A classified systems (col 9: early-type, single lens, clear multiple imaging)
# are kept — see the full table for the excluded L-S-X / E-M-A / E-S-X / E-S-M rows.
_raw = [
    "002927.38+254401.7",
    "020121.39+322829.6",
    "023740.63-064112.9",
    "074249.68+334148.9",
    "075523.52+344539.5",
    "085621.59+201040.5",
    "091859.21+510452.5",
    "111027.11+280838.4",
    "111040.42+364924.4",
    "111634.55+091503.0",
    "114154.71+221628.8",
    "120159.02+474323.2",
    "122656.45+545739.0",
    "222825.76+120503.9",
    "234248.68-012032.5",
]

_COORD_RE = re.compile(r'^(\d{2})(\d{2})(\d{2}\.\d+)([+-])(\d{2})(\d{2})(\d{2}\.\d+)$')


def _build():
    result = {}
    for raw in _raw:
        m = _COORD_RE.match(raw)
        if not m:
            raise ValueError(f"Could not parse coordinate: {raw!r}")
        ra_hh, ra_mm, ra_ss, sign, dec_dd, dec_mm, dec_ss = m.groups()
        name = f"J{ra_hh}{ra_mm}{sign}{dec_dd}{dec_mm}"
        if name in result:
            raise ValueError(f"Duplicate name generated: {name!r} (from {raw!r})")
        result[name] = [f"{ra_hh}:{ra_mm}:{ra_ss}", f"{sign}{dec_dd}:{dec_mm}:{dec_ss}"]
    return result


# Dict mapping lens name -> [RA, Dec] in HH:MM:SS format
# Name format: J{HHMM}_{DDSS} where HHMM = RA hours+minutes, DDSS = Dec degrees+arcseconds
gallery_coords = _build()
