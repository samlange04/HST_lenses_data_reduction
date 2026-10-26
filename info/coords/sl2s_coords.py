import re

# SL2S lens positions, keyed by their SL2S designation (SL2SJhhmmss+ddmmss).
# The designations come from Sonnenfeld et al. 2013 (doi:10.1088/0004-637X/777/2/97), Table 6,
# and, for J0205-0930, Sonnenfeld et al. 2015 (doi:10.1088/0004-637X/800/2/94), Table 3
# (which also named J0233-0438, removed from the roster 2026-10-14). They are TRUNCATED, so the RA part is only good to 1 s of time (up to ~15"): using
# them as positions put J0959+0206's stamp search on a neighbour 2.6" away (2026-10-04 pilot).
# The positions below are the decimal-degree RA/Dec of Gavazzi et al. 2014
# (doi:10.1088/0004-637X/785/2/144), Tables 1 and 3, matched on the FULL designation (the short
# J-name is not unique in that table: J0226-0406 also names 022611-040646). They agree with the
# GO-17130 / GO-11588 target positions to <0.7" (1.5" for one 11588 visit). J0858-0143 is not in
# Gavazzi+14 and is the one SDSS DR16 galaxy inside its designation strip.
# _build() checks that every position lies inside its own designation's truncation strip.
# One table for all four SL2S samples; the roster is in AGENTS.md (*SL2S samples*).
_raw = [
    # sl2s_dinos_a
    ("022610-042011", 36.5444, -4.3366),
    ("023251-040823", 38.2149, -4.1399),
    ("090407-005952", 136.0330, -0.9980),
    ("135949+553550", 209.9567, 55.5973),
    ("140650+522619", 211.7097, 52.4386),
    ("140546+524311", 211.4426, 52.7198),
    ("095921+020638", 149.8407, 2.1107),
    ("021411-040502", 33.5467, -4.0841),
    ("021737-051329", 34.4049, -5.2248),
    # sl2s_dinos_x
    ("020524-093023", 31.3527, -9.5065),
    ("021247-055552", 33.1993, -5.9312),
    ("021801-080247", 34.5053, -8.0465),
    ("022046-094927", 35.1919, -9.8244),
    ("084847-035103", 132.1968, -3.8511),   # J0848-0351, removed 2026-10-24 (no counter-image)
    ("085540-014730", 133.9172, -1.7918),
    ("140454+520024", 211.2269, 52.0068),
    ("220506+014703", 331.2788, 1.7844),
    ("222148+011542", 335.4534, 1.2619),
    ("023307-043838", 38.2794, -4.6440),    # J0233-0438, removed 2026-10-14 (both visits failed)
    ("021902-082934", 34.7589, -8.4930),
    ("084909-041226", 132.2896, -4.2074),
    ("220329+020518", 330.8709, 2.0886),
    # sl2s_wf
    ("021325-074355", 33.3522, -7.7319),    # J0213-0743, removed 2026-10-14 (not a single deflector)
    ("022511-045433", 36.2960, -4.9093),
    ("084959-025142", 132.4990, -2.8619),
    ("085826-014300", 134.610863, -1.716815),  # SDSS DR16 (VizieR V/154): the only galaxy in the designation strip; not in G14
    ("090106-025906", 135.2752, -2.9852),
    ("022648-040610", 36.7016, -4.1029),
    ("142731+551645", 216.8803, 55.2792),
    ("221326-000946", 333.3591, -0.1629),
    ("135847+545913", 209.6958, 54.9870),
    ("140123+555705", 210.3466, 55.9515),
    ("141917+511729", 214.8219, 51.2913),
    # sl2s_snap
    ("140156+554446", 210.4849, 55.7463),
    ("022357-065142", 35.9914, -6.8619),
    ("140221+550534", 210.5897, 55.0930),
    ("141137+565119", 212.9043, 56.8553),
    ("142059+563007", 215.2494, 56.5021),
    ("142321+572243", 215.8413, 57.3786),   # J1423+5722, removed 2026-10-24 (no multiple imaging)
]

_DESIG_RE = re.compile(r'^(\d{2})(\d{2})(\d{2})([+-])(\d{2})(\d{2})(\d{2})$')


def _sexagesimal(deg, hours):
    """'HH:MM:SS.sss' (hours=True) or '+DD:MM:SS.ss' from decimal degrees."""
    v = deg / 15.0 if hours else abs(deg)
    total = round(v * 3600.0, 3 if hours else 2)
    h, rem = divmod(total, 3600.0)
    m, s = divmod(rem, 60.0)
    if hours:
        return f"{int(h):02d}:{int(m):02d}:{s:06.3f}"
    return f"{'-' if deg < 0 else '+'}{int(h):02d}:{int(m):02d}:{s:05.2f}"


def _build():
    result = {}
    for desig, ra, dec in _raw:
        m = _DESIG_RE.match(desig)
        if not m:
            raise ValueError(f"Could not parse designation: {desig!r}")
        hh, mm, ss, sign, dd, dm, ds = m.groups()
        # Truncation strip: RA in [hhmmss, hhmmss + 1 s), |Dec| in [ddmmss, ddmmss + 1"),
        # with a 0.5 s / 1" allowance for catalogue-to-catalogue scatter.
        ra_s = ra / 15.0 * 3600.0 - (int(hh) * 3600 + int(mm) * 60 + int(ss))
        dec_s = abs(dec) * 3600.0 - (int(dd) * 3600 + int(dm) * 60 + int(ds))
        if not (-0.5 <= ra_s < 1.5 and -1.0 <= dec_s < 2.0) or (sign == '-') != (dec < 0):
            raise ValueError(f"{desig}: position ({ra}, {dec}) lies outside its designation "
                             f"strip (RA {ra_s:+.2f} s, |Dec| {dec_s:+.2f} arcsec)")
        name = f"J{hh}{mm}{sign}{dd}{dm}"
        if name in result:
            raise ValueError(f"Duplicate name generated: {name!r} (from {desig!r})")
        result[name] = [_sexagesimal(ra, True), _sexagesimal(dec, False)]
    return result


# Dict mapping lens name -> [RA, Dec] in HH:MM:SS format
sl2s_coords = _build()
