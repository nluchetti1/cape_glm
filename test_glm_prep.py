"""Offline test of glm_prep.py: fake GOES-19 S3 listing serving real GLM-format netCDF files.
    python test_glm_prep.py"""
import datetime as dt, os, sys, tempfile
import numpy as np, netCDF4
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import glm_prep as G
FAILS = []
def check(c, m):
    print(("  ok    " if c else "  FAIL  ") + m)
    if not c: FAILS.append(m)
DAY = dt.datetime(2025, 7, 1)
LC39A = (28.6084, -80.6043)
def nc_bytes(start, pts):
    ds = netCDF4.Dataset("glm.nc", "w", memory=4096, format="NETCDF4")
    ds.time_coverage_start = start.strftime("%Y-%m-%dT%H:%M:%S.0Z")
    ds.createDimension("number_of_flashes", None)
    v = {k: ds.createVariable(k, t, ("number_of_flashes",)) for k, t in
         (("flash_lat", "f4"), ("flash_lon", "f4"), ("flash_time_offset_of_first_event", "f4"), ("flash_quality_flag", "i2"))}
    if pts:
        v["flash_lat"][:] = [p[0] for p in pts]; v["flash_lon"][:] = [p[1] for p in pts]
        v["flash_time_offset_of_first_event"][:] = [p[2] for p in pts]; v["flash_quality_flag"][:] = [p[3] for p in pts]
    return bytes(ds.close())
FILES, KEYS = {}, {}
for h in range(24):
    hr = DAY.replace(hour=h); KEYS[h] = []
    for k in range(180):
        st = hr + dt.timedelta(seconds=20 * k)
        if h == 17 and st.minute >= 30 and k % 2:
            continue                                            # half of 17:30-18:00 missing
        key = f"GLM-L2-LCFA/{st:%Y}/{st.timetuple().tm_yday:03d}/{st:%H}/OR_GLM-L2-LCFA_G19_s{st:%Y}{st.timetuple().tm_yday:03d}{st:%H%M%S}0_e_c.nc"
        pts = []
        if DAY.replace(hour=14, minute=30) <= st < DAY.replace(hour=15):
            pts = [(LC39A[0] + 0.05, LC39A[1], 2.0, 0), (LC39A[0] + 0.05, LC39A[1], 9.0, 0), (LC39A[0] + 0.05, LC39A[1], 15.0, 0),
                   (LC39A[0], LC39A[1], 5.0, 1),             # bad quality: dropped
                   (35.0, -75.0, 5.0, 0)]                     # far outside the map: dropped
        if h == 20 and k == 7:
            FILES[key] = b"not a netcdf file"                  # unreadable: file_ok False
        else:
            FILES[key] = nc_bytes(st, pts)
        KEYS[h].append(key)
class Sess:
    def get(self, url, timeout=None):
        class R:
            def __init__(s, c=b"", t=""): s.content, s.text = c, t
            def raise_for_status(s): pass
        if "list-type=2" in url:
            h = int(url.split("prefix=")[1].split("/")[3])
            return R(t="<ListBucketResult>" + "".join(f"<Contents><Key>{k}</Key></Contents>" for k in KEYS[h]) + "</ListBucketResult>")
        return R(c=FILES[url.split("/", 3)[3]])
G.session = lambda: Sess()
out = tempfile.mkdtemp()
sys.argv = ["glm_prep.py", "2025-07", "--days", "1", "--procs", "1", "--out", out]
G.main()
z = np.load(os.path.join(out, "glm_202507.npz"))
e = lambda t: (t - dt.datetime(1970, 1, 1)).total_seconds()
check(len(z["file_t"]) == 24 * 180 - 45, f"every listed file recorded: {len(z['file_t'])} (24 x 180 minus 45 missing)")
check(int((~z["file_ok"]).sum()) == 1, "the one unreadable file is flagged, not dropped")
check(len(z["flash_t"]) == 90 * 3, f"good flashes in the map kept: {len(z['flash_t'])} (bad-quality and far-away dropped)")
check(((z["flash_t"] > e(DAY.replace(hour=14, minute=30))) & (z["flash_t"] <= e(DAY.replace(hour=15)))).all(), "flash times are file start + offset")
check(set(z["satellite"].tolist()) == {"G19"} and z["flash_lat"].dtype == np.float32, "satellite and dtypes recorded")
import shutil; shutil.copy(os.path.join(out, "glm_202507.npz"), "/tmp/glm_202507.npz")
print(f"\n{len(FAILS)} failure(s)"); sys.exit(1 if FAILS else 0)
