#!/usr/bin/env python3
"""
Cape GLM archive - extract one month of GOES-East GLM flashes over the CloudScope map, once.

The climatologies need "any GLM flash within 10 / 20 nmi of a pad in the 30 minutes before each
hourly sample". Reading the ~90 twenty-second GLM files per sample inside every climatology run cost
~2 h per month job (Oct 2026: a month ran past 3 h). This job reads every file of a month ONCE and
keeps only what the climatologies use:

    flash_t    float64  flash time, seconds since 1970-01-01 (first event)
    flash_lat  float32
    flash_lon  float32
    file_t     int64    start time of every GLM file the archive lists, seconds since 1970
    file_ok    bool     whether that file was downloaded and parsed
    satellite  str      per-file "G16" / "G19"

Same filter as the live classifier's reader (mrms_classify._glm_read): flash_quality_flag == 0, and
inside the map (27.6-29.6 N, 81.6-79.6 W) widened by 25 nmi. GOES-16 before 2025-04-07, GOES-19
after (the climatologies' switch date). Files are downloaded on threads and parsed serially:
netCDF4/HDF5 is not thread-safe (parsing on download threads segfaulted in testing).

    python glm_prep.py 2021-07 --procs 4 --out out        -> out/glm_202107.npz
"""
import argparse
import calendar
import concurrent.futures as cf
import datetime as dt
import logging
import multiprocessing as mp
import os
import re
import tempfile
import time

import numpy as np
import requests

DOMAIN = {"lat_min": 27.6, "lat_max": 29.6, "lon_min": -81.6, "lon_max": -79.6}   # = mrms_classify
MARGIN_NM = 25.0                                                                     # = GLM_MARGIN_NM
GOES_SWITCH = dt.datetime(2025, 4, 7)
BUCKETS = ("https://noaa-goes16.s3.amazonaws.com", "https://noaa-goes19.s3.amazonaws.com")
THREADS = 16
EPOCH = dt.datetime(1970, 1, 1)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(processName)s %(message)s")


def bucket(t):
    return BUCKETS[0] if t < GOES_SWITCH else BUCKETS[1]


def session():
    s = requests.Session()
    s.mount("https://", requests.adapters.HTTPAdapter(max_retries=3, pool_maxsize=THREADS * 2))
    s.headers["User-Agent"] = "CloudScope-GLM-archive/1.0"
    return s


def hour_keys(sess, hr):
    """[(file start, key)] in one GLM-L2-LCFA hour folder."""
    out, token = [], None
    prefix = f"GLM-L2-LCFA/{hr:%Y}/{hr.timetuple().tm_yday:03d}/{hr:%H}/"
    while True:
        url = (f"{bucket(hr)}/?list-type=2&prefix={prefix}"
               + (f"&continuation-token={requests.utils.quote(token)}" if token else ""))
        r = sess.get(url, timeout=60)
        r.raise_for_status()
        for key in re.findall(r"<Key>([^<]+)</Key>", r.text):
            m = re.search(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})", key)
            if m:
                y, d, H, M, S = map(int, m.groups())
                out.append((dt.datetime(y, 1, 1, H, M, S) + dt.timedelta(days=d - 1), key))
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", r.text)
        if not m:
            return sorted(set(out))
        token = m.group(1)


def read_file(raw):
    """(lat, lon, epoch seconds) of good flashes inside the map + margin - _glm_read's filter."""
    import netCDF4
    try:
        ds = netCDF4.Dataset("glm", memory=raw)
    except Exception:
        fd, path = tempfile.mkstemp(suffix=".nc")
        with os.fdopen(fd, "wb") as f:
            f.write(raw)
        ds = netCDF4.Dataset(path)
        os.remove(path)
    try:
        lat = np.asarray(ds.variables["flash_lat"][:], float)
        lon = np.asarray(ds.variables["flash_lon"][:], float)
        off = np.asarray(ds.variables["flash_time_offset_of_first_event"][:], float)
        q = (np.asarray(ds.variables["flash_quality_flag"][:], int)
             if "flash_quality_flag" in ds.variables else np.zeros(lat.shape, int))
        base = dt.datetime.strptime(ds.time_coverage_start[:19], "%Y-%m-%dT%H:%M:%S")
    finally:
        ds.close()
    mlat, mlon = MARGIN_NM * 1.852 / 111.32, MARGIN_NM * 1.852 / (111.32 * 0.87)
    keep = ((q == 0) & (lat >= DOMAIN["lat_min"] - mlat) & (lat <= DOMAIN["lat_max"] + mlat)
            & (lon >= DOMAIN["lon_min"] - mlon) & (lon <= DOMAIN["lon_max"] + mlon))
    t0 = (base - EPOCH).total_seconds()
    return lat[keep], lon[keep], t0 + off[keep]


def fetch(sess, url):
    try:
        r = sess.get(url, timeout=60)
        r.raise_for_status()
        return r.content
    except Exception:
        return None


def do_hours(hours):
    """Every GLM file in these hour folders. Returns arrays (see module docstring)."""
    sess = session()
    FT, FLA, FLO, KT, KO, KS = [], [], [], [], [], []
    t_start, nfiles, tdl, tparse = time.monotonic(), 0, 0.0, 0.0
    pool = cf.ThreadPoolExecutor(max_workers=THREADS)
    for n, hr in enumerate(hours):
        try:
            keys = hour_keys(sess, hr)
        except Exception as e:
            logging.warning(f"{hr:%Y-%m-%d %HZ} listing: {type(e).__name__}: {e}")
            keys = []
        a = time.monotonic()
        raws = list(pool.map(lambda k: fetch(sess, f"{bucket(hr)}/{k[1]}"), keys))
        b = time.monotonic()
        for (st, key), raw in zip(keys, raws):
            ok = False
            if raw is not None:
                try:
                    la, lo, tt = read_file(raw)
                    FLA.append(la.astype(np.float32)); FLO.append(lo.astype(np.float32)); FT.append(tt)
                    ok = True
                except Exception:
                    pass
            KT.append(int((st - EPOCH).total_seconds())); KO.append(ok)
            KS.append("G16" if "_G16_" in key else "G19" if "_G19_" in key else "?")
        tdl += b - a
        tparse += time.monotonic() - b
        nfiles += len(keys)
        if (n + 1) % 24 == 0 or n + 1 == len(hours):
            el = time.monotonic() - t_start
            logging.info(f"{hr:%Y-%m-%d %HZ}: {n + 1}/{len(hours)} hours, {nfiles} files, "
                         f"{nfiles / max(el, 1e-9):.1f} files/s (download {tdl:.0f} s, parse {tparse:.0f} s)")
    pool.shutdown()
    cat = lambda L, d: np.concatenate(L).astype(d) if L else np.zeros(0, d)
    return {"flash_t": cat(FT, np.float64), "flash_lat": cat(FLA, np.float32), "flash_lon": cat(FLO, np.float32),
            "file_t": np.array(KT, np.int64), "file_ok": np.array(KO, bool), "satellite": np.array(KS)}


def month_hours(ym, days=None):
    y, m = map(int, ym.split("-"))
    nd = days or calendar.monthrange(y, m)[1]
    return [dt.datetime(y, m, 1) + dt.timedelta(hours=h) for h in range(24 * nd)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("month")
    ap.add_argument("--procs", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--days", type=int, default=None, help="first N days only (testing)")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    hours = month_hours(a.month, a.days)
    k = max(1, min(a.procs, len(hours)))
    chunks = [hours[i * len(hours) // k:(i + 1) * len(hours) // k] for i in range(k)]
    t0 = time.monotonic()
    if len(chunks) == 1:
        parts = [do_hours(chunks[0])]
    else:
        with mp.get_context("spawn").Pool(len(chunks)) as pool:
            parts = pool.map(do_hours, chunks)
    res = {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"glm_{a.month.replace('-', '')}.npz")
    np.savez_compressed(path, **res, month=np.array(a.month), hours=np.array(len(hours)))
    logging.info(f"{a.month}: {len(res['file_t'])} files ({int(res['file_ok'].sum())} read), "
                 f"{len(res['flash_t'])} flashes kept, {time.monotonic() - t0:.0f} s -> {path} "
                 f"({os.path.getsize(path) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
