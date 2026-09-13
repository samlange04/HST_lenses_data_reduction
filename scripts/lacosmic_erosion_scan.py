"""Read-only LACosmic sigclip/objlim scan: is CR rejection eating real signal?

Nothing is written -- FLC DQ is read, never updated, so this is safe to run against
the live calibrated tree while a campaign is in flight.

__The metric, and why the obvious one does not work__

The first thing to try is "pixels flagged as CR in EVERY frame at the same sky
position" -- a real cosmic ray hits one frame, so a position flagged in all of them
must be real signal. It sounds right and it DOES NOT WORK: run it on J1420+6019
F814W, the one product documented to have been eroded at the default 4.5/5.0, and it
returns zero coincidences at every setting. Dither means each frame samples the sky
at a different pixel phase, so the flagged pixels never line up exactly. Do not
resurrect that metric without re-checking it against the control below.

What works instead: project each frame's CR flags onto the NO-CR drizzle (which has
had nothing removed, so it is the truth map) and ask what fraction of the genuinely
bright pixels near the lens get flagged.

__Read the RESPONSE to the thresholds, not the absolute number__

The absolute fraction is not comparable across instruments or lenses -- it scales
with CR rate, frame count and pixel scale, and on gallery the best-behaved lens
scores the highest. Erosion is over-aggressive thresholds, so it RELAXES when they
are loosened; a genuine CR sits far above threshold and is flagged regardless.

    J1420+6019 F814W (ACS, documented erosion)  2.8% -> 1.3%   -54%   <-- erosion
    J2228+1205 / J0237-0641 F606W (gallery)      5.3/14.9% -> 4.6/13.4%   -13/-10%
    J1201+4743 / J1110+3649 F606W (gallery)     15.9/11.8% -> 11.5/7.7%   -28/-35%

The ACS control halves; the gallery lenses barely move, and the flagged pair moves
LESS than the well-behaved pair. That is what "no tuning would help" looks like.

Usage:
    uv run python scripts/lacosmic_erosion_scan.py <sample> <lens> <filt> <acs|uvis>
"""
import sys, warnings, glob, numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, 'info')
import astroscrappy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import astropy.units as u

def run(sample, lens, filt, instr, combos):
    if instr == 'acs':
        rd, sat, pre = 4.5, 84700.0, 'acs_wfc_flc'
    else:
        rd, sat, pre = 3.1, 70000.0, 'wfc3_uvis_flc'
    if sample == 'gallery':
        from gallery_coords import gallery_coords as cat
    else:
        from slacs_coords import slacs_coords as cat
    c = SkyCoord(*cat[lens], unit=(u.hourangle, u.deg))
    nc = glob.glob(f'data/drizzled/{sample}/{lens}/{filt}/*_nocrrej_drc_sci.fits')[0]
    D = fits.getdata(nc).astype(float); H = fits.getheader(nc); W = WCS(H)
    _, md, sd = sigma_clipped_stats(D[np.isfinite(D)])
    SN = (D - md) / sd
    px, py = W.world_to_pixel(c)
    ps = abs(H.get('CD1_1', 1e-5)) * 3600
    Y, X = np.mgrid[:D.shape[0], :D.shape[1]]
    near = np.hypot(Y - py, X - px) * ps < 2.0
    truth = near & (SN > 10)                       # real, bright, at the lens
    print(f'=== {lens} {filt} ({instr}) -- {int(truth.sum())} high-S/N sky pixels within 2" ===')
    files = sorted(glob.glob(f'data/calibrated/{sample}/{lens}/{filt}/*_flc.fits'))
    frames = []
    for f in files:
        with fits.open(f) as h:
            for ext in (1, 2):
                frames.append((h['SCI', ext].data.astype(np.float32),
                               (h['DQ', ext].data & (4 | 8 | 128 | 512)) > 0,
                               WCS(h['SCI', ext].header, h)))
    for sig, obj in combos:
        hit = np.zeros(D.shape, bool); tot = 0
        for sci, bad, w in frames:
            m, _ = astroscrappy.detect_cosmics(
                sci, inmask=bad, sigclip=sig, sigfrac=0.3, objlim=obj, gain=1.0,
                readnoise=rd, satlevel=sat, niter=4, sepmed=True,
                cleantype='medmask', fsmode='median')
            yy, xx = np.nonzero(m)
            if not len(yy): continue
            sk = w.pixel_to_world(xx, yy)
            iy, ix = W.world_to_pixel(sk)
            iy = np.round(ix).astype(int), np.round(iy).astype(int)
            ok = ((iy[0] >= 0) & (iy[0] < D.shape[0]) & (iy[1] >= 0) & (iy[1] < D.shape[1]))
            hit[iy[0][ok], iy[1][ok]] = True
            tot += int(ok.sum())
        eroded = int((hit & truth).sum())
        frac = eroded / max(1, int(truth.sum()))
        print(f'  sigclip {sig:>4} objlim {obj:>4} | CR flags total {tot:>6} | '
              f'landing on real signal: {eroded:>4} = {frac:>6.1%} of the source')

if __name__ == '__main__':
    run(*sys.argv[1:5], [(4.5,5.0),(6.0,5.0),(8.0,8.0),(10.0,12.0)])
