"""Drop GSC240-only visits when another visit of the same lens+band has a fitted WCS.

A visit whose frames carry only a GSC240 guide-star WCS (no FIT_* astrometric fit) does not
agree in absolute position with a visit that has a fitted Gaia/GSC242/HSC solution, so
combining them smears and shifts the deflector:
  * J0959+0206 F475X (WFC3/UVIS, GO-17130): the GSC240 visit sat ~0.5" from the Gaia-fit
    visit, which itself agreed with ACS F814W to ~0.05"; combined, the deflector moved 0.14".
  * J0912+0029 F814W (ACS, GO-10494): visit 02 (2005, GSC240) and visit 52 (2006,
    FIT_REL_GAIAeDR3). The tracked stamp is visit 52 alone; the 2026-10-05 re-drizzle combined
    both because the ACS driver lacked this rule.
A lens whose every visit is GSC240-only keeps them all (there is nothing better to choose).
Used by drizzle_wfc3_uvis.py and drizzle_acs_wfc.py, which share the <rootname>_flc.fits layout.
"""
import os

from astropy.io import fits


def visit_wcs_drop(files):
    """Rootnames to exclude: frames of GSC240-only visits when a fitted visit exists.
    A visit is the 6-character rootname prefix (ipppss)."""
    by_visit = {}
    for f in files:
        root = os.path.basename(f).split('_')[0]
        name = str(fits.getval(f, 'WCSNAME', ext=('SCI', 1)))
        by_visit.setdefault(root[:6], []).append((root, 'FIT' in name.upper()))
    fitted = {v for v, fr in by_visit.items() if all(ok for _, ok in fr)}
    if not fitted:
        return set()
    return {root for v, fr in by_visit.items() if v not in fitted for root, _ in fr}
