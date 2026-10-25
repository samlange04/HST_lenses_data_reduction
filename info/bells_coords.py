import re

# Source: Brownstein et al. 2012 (doi:10.1088/0004-637X/744/1/41), Table 2, column 1
# ("BOSS Properties of BELLS Candidate Systems")
# Raw format: HHMMSS.ss[+/-]DDMMSS.s (J2000)
# Only grade-A single-deflector systems are kept: Table 3 classification E-S-A or L-S-A
# (21 + 3), the same rule as slacs_other. Excluded grade-A row: J1522+2910 (E-M-A).
# Three of the 24 were later removed from the roster (info/lens_samples.json) and are
# marked below; their positions stay so the designation list matches the paper.
_raw = [
    "015107.37+004909.0",
    "074734.75+444859.3",
    "074724.12+505537.5",
    "080105.30+472749.6",
    "083049.73+511631.8",
    "094427.47-014742.4",
    "115944.63-000728.2",
    "121504.44+004726.0",
    "122151.92+380610.5",
    "123427.99-024129.6",
    "131829.39-010421.6",
    "133751.31+362018.1",   # J1337+3620, removed 2026-10-14 (companion needs its own mass)
    "134910.30+361239.7",   # J1349+3612, removed 2026-10-14 (same-z companion)
    "135218.99+321651.8",
    "154118.56+181235.1",
    "154246.33+162951.8",   # J1542+1629, removed 2026-10-22 (neighbour mask reaches theta_E)
    "154503.57+274805.3",
    "160113.27+213833.9",
    "161109.80+170526.6",
    "163150.33+185404.1",
    "163714.58+143930.1",
    "212252.04+040935.5",
    "212510.67+041131.6",
    "230335.17+003703.2",]

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
# Name format: J{HHMM}{+/-DDMM} where HHMM = RA hours+minutes, DDMM = Dec degrees+arcminutes
bells_coords = _build()
