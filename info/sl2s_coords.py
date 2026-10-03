import re

# Source: the SL2S designations themselves (SL2SJhhmmss+ddmmss, truncated J2000), as listed in
# Sonnenfeld et al. 2013 (doi:10.1088/0004-637X/777/2/97), Table 6, and, for J0205-0930 and
# J0233-0438, Sonnenfeld et al. 2015 (doi:10.1088/0004-637X/800/2/94), Table 3. Neither paper
# tabulates more precise positions, so these are good to ~1": ample for the common-WCS tangent
# point, and make_cutouts.py re-centres on the deflector anyway.
# One table for all four SL2S samples (sl2s_dinos_a, sl2s_dinos_x, sl2s_pc, sl2s_snap); the
# roster and the reasons for it are in AGENTS.md (*SL2S samples*).
_raw = [
    # sl2s_dinos_a
    "022610-042011", "023251-040823", "090407-005952", "135949+553550", "140650+522619",
    "140546+524311", "095921+020638", "021411-040502", "021737-051329",
    # sl2s_dinos_x
    "020524-093023", "021247-055552", "021801-080247", "022046-094927", "084847-035103",
    "085540-014730", "140454+520024", "220506+014703", "222148+011542", "023307-043838",
    "021902-082934", "084909-041226", "220329+020518",
    # sl2s_pc
    "021325-074355", "022511-045433", "084959-025142", "085826-014300", "090106-025906",
    "022648-040610", "142731+551645", "221326-000946", "135847+545913", "140123+555705",
    "140156+554446", "141917+511729",
    # sl2s_snap
    "022357-065142", "140221+550534", "141137+565119", "142059+563007", "142321+572243",
]

_COORD_RE = re.compile(r'^(\d{2})(\d{2})(\d{2})([+-])(\d{2})(\d{2})(\d{2})$')


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
sl2s_coords = _build()
