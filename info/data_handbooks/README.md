# Data handbooks and PSF references

The STScI handbooks and papers that the pipeline's DQ-bit, pixel-scale, read-noise and
PSF choices cite (see `instrument_drizzle_ref.md`, the condensed version with page
numbers). **The PDFs are not tracked** (147 MB, untracked 2026-10-04; `*.pdf` here is
gitignored). Drop them into this directory under the filenames below so the page
citations in `instrument_drizzle_ref.md` and the scripts resolve.

| file | document | source |
|---|---|---|
| `acs_dhb.pdf` | ACS Data Handbook v14.0 (Cycle 34, 2026), 240 pp. | https://hst-docs.stsci.edu/acsdhb |
| `wfc3dhb2024_final.pdf` | WFC3 Data Handbook v6.0 (June 2024), 296 pp. | https://hst-docs.stsci.edu/wfc3dhb |
| `The_DrizzlePac_Handbook_Version3.pdf` | The DrizzlePac Handbook v3.0, 154 pp. | https://hst-docs.stsci.edu/drizzpac |
| `wfpc2_dhb.pdf` | WFPC2 Data Handbook v5.0 (July 2010), 176 pp. | https://www.stsci.edu/hst/instrumentation/legacy/wfpc2 (Documentation) |
| `nic_dhb.pdf` | NICMOS Data Handbook v8.0 (May 2009), 187 pp. | https://www.stsci.edu/hst/instrumentation/legacy/nicmos (Documentation) |
| `wfc3_ir_psfs.pdf` | Anderson 2016, "Empirical Models for the WFC3/IR PSF", WFC3 ISR 2016-12 | https://www.stsci.edu/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs |
| `wfpc2_wfc3_database.pdf` | Dauphin et al. 2021, "The WFPC2 and WFC3 PSF Database", WFC3 ISR 2021-12 | same ISR index |
| `wfpc2_psf_2000.pdf` | Anderson & King 2000, "Toward High-Precision Astrometry with WFPC2. I. Deriving an Accurate PSF", PASP 112, 1360 | https://arxiv.org/abs/astro-ph/0006325 |
| `HSTphot.pdf` | Dolphin 2000, "WFPC2 Stellar Photometry with HSTphot", PASP 112, 1383 | https://arxiv.org/abs/astro-ph/0006217 |

The handbooks are living documents: the versions above are the ones the page numbers in
this repo were read from. A newer edition may renumber pages.
