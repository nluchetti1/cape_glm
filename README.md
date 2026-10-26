# Cape GLM archive

GOES-East GLM flashes over the CloudScope map, extracted ONCE so the climatologies never read the
raw GLM archive. One file per month at `https://<owner>.github.io/cape_glm/data/glm_YYYYMM.npz`:
flash time / lat / lon (quality flag 0, map + 25 nmi - the live classifier's own filter), plus the
start time and read status of every GLM file, so a climatology can tell whether a sample's
30-minute window was complete. GOES-16 before 2025-04-07, GOES-19 after.

## Files
```
cape_glm/
├── glm_prep.py              extracts one month
├── requirements.txt         numpy, requests, netCDF4
├── test_glm_prep.py         offline test (fake S3, real-format netCDF)
└── .github/workflows/glm.yml
```

## Setup and first run
1. Create the repo `cape_glm` (public) and commit these files. The climatologies look for
   `https://nluchetti1.github.io/cape_glm/data` - keep that name, or set `GLM_URL` in climo_run.py.
2. Actions -> "Cape GLM archive" -> Run: start `2020-10`, end `2025-12`, rerun unchecked.
   ~130,000 GLM files per month; the log prints files/s every 24 hours of data.
3. Settings -> Pages -> deploy from `gh-pages`. The index page lists every month with its file
   count, % read and flashes kept; `data/index.json` has the same.
4. Only then switch the climatologies to their archive-reading climo_run.py (after any run in
   progress has finished). No re-score is needed: the flashes and the completeness rule are the same.
