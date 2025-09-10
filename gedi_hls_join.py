import os, sys, io, json, warnings, xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from dateutil import parser as dtparser
from typing import List, Tuple
import tempfile, os
import h5py
import re

import requests
import earthaccess
import netCDF4 as nc
import numpy as np
import pandas as pd
from shapely.geometry import box
from tqdm import tqdm

import rasterio
from rasterio.errors import RasterioIOError
from pystac_client import Client
import planetary_computer as pc

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from dateutil.relativedelta import relativedelta  # you already have this
from urllib.parse import urlparse

import argparse
import json
from shapely.geometry import shape
from shapely import wkt as shapely_wkt

# ---------------- CLI / CONFIG ----------------
def _parse_months(mstr: str) -> list[int]:
    mstr = (mstr or "").strip().lower()
    if mstr in ("all", "1-12"):
        return list(range(1, 13))
    if not mstr:
        return [7]  # sensible default
    return [int(x) for x in mstr.split(",") if x.strip()]

def _resolve_aoi_bbox(args) -> tuple[float,float,float,float]:
    # Priority: --aoi-bbox > --aoi-geojson > --aoi-wkt > default
    if args.aoi_bbox:
        parts = [float(x) for x in args.aoi_bbox.split(",")]
        if len(parts) != 4:
            raise ValueError("--aoi-bbox needs 4 comma-separated numbers: minLon,minLat,maxLon,maxLat")
        return tuple(parts)  # type: ignore[return-value]

    if args.aoi_geojson:
        with open(args.aoi_geojson, "r") as f:
            gj = json.load(f)
        if gj.get("type") == "FeatureCollection":
            geom = shape(gj["features"][0]["geometry"])
        elif gj.get("type") == "Feature":
            geom = shape(gj["geometry"])
        else:
            geom = shape(gj)  # bare geometry dict
        minx, miny, maxx, maxy = geom.bounds
        return (float(minx), float(miny), float(maxx), float(maxy))

    if args.aoi_wkt:
        geom = shapely_wkt.loads(args.aoi_wkt)
        minx, miny, maxx, maxy = geom.bounds
        return (float(minx), float(miny), float(maxx), float(maxy))

    # Default AOI (your current one)
    return (-97.5, 46.5, -90.0, 49.5)

def parse_cli():
    p = argparse.ArgumentParser(
        description="Join GEDI L4A to HLS (L30+S30) at GEDI shot locations and export CSV."
    )
    # AOI
    p.add_argument("--aoi-bbox", type=str, help="minLon,minLat,maxLon,maxLat")
    p.add_argument("--aoi-geojson", type=str, help="Path to GeoJSON (Feature/FeatureCollection/Geometry)")
    p.add_argument("--aoi-wkt", type=str, help="WKT geometry string")

    # Time selection
    p.add_argument("--year", type=int, default=2022, help="Target year (used with --months)")
    p.add_argument("--months", type=str, default="7",
                   help="Comma-separated months (e.g. '6,7,8') or 'all'")

    # Rolling window & limits
    p.add_argument("--rolling-days", type=int, default=7,
                   help="+/- days around each month for HLS median")
    p.add_argument("--max-hls-per-month", type=int, default=200,
                   help="Max HLS granules (L30+S30 combined) per month")
    p.add_argument("--limit-gedi-per-month", type=int, default=40,
                   help="Limit GEDI granules per month (None to read all)")

    # Output & workspace
    p.add_argument("--workdir", type=str, default=".",
                   help="Folder where caches & outputs are stored")
    p.add_argument("--out-csv", type=str, default="gedi_hls_monthly_30m.csv",
                   help="Output CSV (relative to --workdir if not absolute)")

    # Product/quality options
    p.add_argument("--gedi-collection-id", default="C2237824918-ORNL_CLOUD")
    p.add_argument("--accept-water", action="store_true",
                   help="Keep water pixels flagged as clear water")
    p.add_argument("--accept-adjacent", action="store_true",
                   help="Keep pixels adjacent to cloud/shadow")
    p.add_argument("--accept-snow", action="store_true",
                   help="Keep snow/ice pixels")

    return p.parse_args()


# -------------- Helpers --------------
def looks_like_netcdf(buf: bytes) -> bool:
    # NetCDF-3 starts with b"CDF", NetCDF-4 is HDF5 and starts with b"\x89HDF\r\n\x1a\n"
    return (len(buf) >= 3 and buf[:3] == b"CDF") or (len(buf) >= 8 and buf[:8] == b"\x89HDF\r\n\x1a\n")

def _granule_name_with_prefix(fname: str, short_name: str) -> str:
    # ensure .h5 suffix
    if not fname.lower().endswith(".h5"):
        fname = fname + ".h5"
    # ensure product prefix
    if not fname.startswith(short_name + "."):
        fname = f"{short_name}.{fname}"
    return fname

def resolve_hyrax_base_url(sess: requests.Session, base_url: str,
                           collection_id: str, short_name: str,
                           timeout: int = 30) -> str | None:
    """
    Given a guessed Hyrax 'base_url', try URL variants until .dmr.xml returns 200.
    Returns the working base URL (no suffix) or None if none worked.
    """
    u = urlparse(base_url)
    # extract filename
    fname = u.path.rsplit("/", 1)[-1]  # e.g., "GEDI04_A_....h5" or already with prefix
    fname_pref = _granule_name_with_prefix(fname, short_name)

    # candidate base URLs to try (no suffix)
    cands = []

    # 1) what you already have
    cands.append(base_url.rstrip(".dmr").rstrip(".dmr.xml").rstrip(".dmr.html"))

    # 2) same path but replace filename with prefixed variant
    cands.append(base_url.replace(fname, fname_pref))

    # 3) add providers segment (no prefix)
    cands.append(f"https://opendap.earthdata.nasa.gov/providers/ORNL_CLOUD/collections/{collection_id}/granules/{fname}")

    # 4) add providers + prefixed filename
    cands.append(f"https://opendap.earthdata.nasa.gov/providers/ORNL_CLOUD/collections/{collection_id}/granules/{fname_pref}")

    seen = set()
    for b in cands:
        if b in seen:
            continue
        seen.add(b)
        try:
            r = sess.get(b + ".dmr.xml", timeout=timeout)
            if r.status_code == 200:
                return b  # success!
        except Exception:
            pass
    return None


def month_timerange(year: int, month: int, rolling_days: int):
    m_start = datetime(year, month, 1)
    m_end = (m_start + relativedelta(months=1)) - timedelta(days=1)
    q_start = m_start - timedelta(days=rolling_days)
    q_end   = m_end + timedelta(days=rolling_days)
    return m_start, m_end, q_start, q_end

def cmr_search_granule_opendap_urls(collection_id: str,
                                    bbox: Tuple[float,float,float,float],
                                    start: datetime, end: datetime,
                                    page_size: int = 2000) -> List[str]:
    """
    Search CMR (UMM-JSON) and extract Hyrax OPeNDAP URLs from RelatedUrls.
    """
    base = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
    params = {
        "collection_concept_id": collection_id,
        "temporal": f"{start.strftime('%Y-%m-%dT%H:%M:%SZ')},{end.strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "bounding_box": f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}",
        "page_size": page_size,
        "page_num": 1
    }
    urls = []
    while True:
        r = requests.get(base, params=params, timeout=60)
        r.raise_for_status()
        js = r.json()
        items = js.get("items", [])
        if not items: break
        for it in items:
            umm = it.get("umm", {})
            opendap_url = None
            for ru in umm.get("RelatedUrls", []):
                u = ru.get("URL")
                if u and "opendap.earthdata.nasa.gov" in u:
                    opendap_url = u.replace(".dmr.html","").replace(".dmr","")
                    break
            if opendap_url:
                urls.append(opendap_url)
        if len(items) < page_size: break
        params["page_num"] += 1
    # dedupe
    seen=set(); uniq=[]
    for u in urls:
        if u not in seen:
            uniq.append(u); seen.add(u)
    return uniq

def cmr_search_granule_opendap_urls_chunked(collection_id: str,
                                            bbox: tuple[float,float,float,float],
                                            start: datetime, end: datetime,
                                            page_size: int = 200,
                                            read_timeout: int = 180) -> list[str]:
    """
    Robust CMR UMM-JSON search:
      - splits the temporal window into 1-month chunks
      - uses HTTP retries with backoff
      - extracts Hyrax OPeNDAP URLs from RelatedUrls
    """
    base = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"

    sess = requests.Session()
    retry = Retry(
        total=5, connect=5, read=5, status=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET"}
    )
    sess.mount("https://", HTTPAdapter(max_retries=retry))

    def temporal_str(a: datetime, b: datetime):
        return f"{a.strftime('%Y-%m-%dT%H:%M:%SZ')},{b.strftime('%Y-%m-%dT%H:%M:%SZ')}"

    urls: list[str] = []
    cur = start
    while cur <= end:
        chunk_end = (cur + relativedelta(months=1)) - timedelta(seconds=1)
        if chunk_end > end:
            chunk_end = end

        page = 1
        while True:
            params = {
                "collection_concept_id": collection_id,
                "temporal": temporal_str(cur, chunk_end),
                "bounding_box": f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}",
                "page_size": page_size,
                "page_num": page,
            }
            try:
                r = sess.get(base, params=params, timeout=read_timeout)
                r.raise_for_status()
            except Exception as e:
                warnings.warn(f"CMR {cur:%Y-%m} page {page} failed: {e}")
                break

            items = r.json().get("items", [])
            if not items:
                break

            for it in items:
                umm = it.get("umm", {})
                opendap_url = None
                for ru in umm.get("RelatedUrls", []):
                    u = ru.get("URL")
                    if u and "opendap.earthdata.nasa.gov" in u:
                        # normalize so netCDF/bytes requests work
                        opendap_url = u.replace(".dmr.html", "").replace(".dmr", "")
                        break
                if opendap_url:
                    urls.append(opendap_url)

            if len(items) < page_size:
                break
            page += 1

        cur = chunk_end + timedelta(seconds=1)

    # dedupe, keep order
    seen = set(); uniq = []
    for u in urls:
        if u not in seen:
            uniq.append(u); seen.add(u)
    return uniq

def ea_search_granule_opendap_urls_chunked(collection_id, bbox, start, end):
    """
    Use earthaccess to query CMR month-by-month and extract Hyrax OPeNDAP URLs.
    If a granule doesn't publish an explicit OPeNDAP RelatedUrl, fall back to
    constructing the Hyrax URL from its GranuleUR.
    """
    urls = []
    cur = start
    while cur <= end:
        chunk_end = (cur + relativedelta(months=1)) - timedelta(seconds=1)
        if chunk_end > end:
            chunk_end = end

        # Query month with concept_id, fallback to short_name if concept_id fails
        try:
            results = earthaccess.search_data(
                concept_id=collection_id,
                bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),  # tuple!
                temporal=(cur, chunk_end)
            )
        except Exception:
            results = earthaccess.search_data(
                short_name="GEDI_L4A_AGB_Density_V2_1",
                bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),
                temporal=(cur, chunk_end)
            )

        found = 0
        for g in results:
            umm = getattr(g, "umm", {}) or {}
            opendap_url = None

            # 1) Try explicit RelatedUrls OPeNDAP link
            for ru in umm.get("RelatedUrls", []) or []:
                u = ru.get("URL")
                if not u:
                    continue
                if "opendap.earthdata.nasa.gov" in u:
                    opendap_url = u.replace(".dmr.html", "").replace(".dmr", "")
                    break

            # 2) Fallback: build from GranuleUR
            if not opendap_url:
                granule_ur = umm.get("GranuleUR")
                if granule_ur:
                    if not granule_ur.endswith(".h5"):
                        granule_ur = granule_ur + ".h5"
                    opendap_url = (
                        f"https://opendap.earthdata.nasa.gov/collections/"
                        f"{collection_id}/granules/{granule_ur}"
                    )

            if opendap_url:
                urls.append(opendap_url); found += 1

        print(f"  {cur:%Y-%m}: granules={len(results)}, opendap_urls+={found}, total={len(urls)}")
        cur = chunk_end + timedelta(seconds=1)

    # dedupe, keep order
    seen = set(); uniq = []
    for u in urls:
        if u not in seen:
            uniq.append(u); seen.add(u)
    return uniq

def ea_search_granule_opendap_urls_for_months(collection_id, bbox, year, months, rolling_days, max_per_month=None):
    """
    Search only the requested months. For each granule:
      1) Use OPeNDAP link in RelatedUrls if present.
      2) Else use UMM['GranuleUR'] if present.
      3) Else derive the filename from a regular HTTPS data link (granule.data_links()).
    Returns a de-duplicated list of Hyrax base URLs (no .dmr/.dmr.html suffix).
    """
    urls = []
    for m in sorted(set(months)):
        m_start = datetime(year, m, 1)
        m_end   = (m_start + relativedelta(months=1)) - timedelta(days=1)
        q_start = m_start - timedelta(days=rolling_days)
        q_end   = m_end + timedelta(days=rolling_days)

        # concept_id first; fallback to short_name
        try:
            results = earthaccess.search_data(
                concept_id=collection_id,
                bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),
                temporal=(q_start, q_end)
            )
        except Exception:
            results = earthaccess.search_data(
                short_name="GEDI_L4A_AGB_Density_V2_1",
                bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),
                temporal=(q_start, q_end)
            )

        if max_per_month:
            results = results[:max_per_month]

        found = 0
        for g in results:
            umm = getattr(g, "umm", {}) or {}
            opendap_url = None

            # (1) explicit OPeNDAP RelatedUrl
            for ru in umm.get("RelatedUrls", []) or []:
                u = ru.get("URL")
                if u and "opendap.earthdata.nasa.gov" in u:
                    opendap_url = u.replace(".dmr.html", "").replace(".dmr", "")
                    break

            # (2) GranuleUR fallback
            if not opendap_url:
                granule_ur = umm.get("GranuleUR")

                # (3) Filename from HTTPS data link if GranuleUR missing
                if not granule_ur:
                    try:
                        for dl in g.data_links() or []:
                            # prefer regular https (not s3) and .h5 files
                            if dl.startswith("http") and dl.lower().endswith(".h5"):
                                granule_ur = dl.rsplit("/", 1)[-1]
                                break
                    except Exception:
                        pass

                if granule_ur:
                    if not granule_ur.endswith(".h5"):
                        granule_ur += ".h5"
                    opendap_url = (
                        f"https://opendap.earthdata.nasa.gov/collections/"
                        f"{collection_id}/granules/{granule_ur}"
                    )

            if opendap_url:
                urls.append(opendap_url); found += 1

        print(f"  {year}-{m:02d}: granules={len(results)}, opendap_urls+={found}, total={len(urls)}")

    # de-dup, keep order
    seen, uniq = set(), []
    for u in urls:
        if u not in seen:
            uniq.append(u); seen.add(u)
    return uniq


def list_beams_from_dmr_xml(sess: requests.Session, base_url: str) -> List[str]:
    """
    Fetch DMR XML and return group names that start with 'BEAM'.
    """
    r = sess.get(base_url + ".dmr.xml", timeout=60)
    r.raise_for_status()
    root = ET.fromstring(r.text)
    beams = set()
    # Hyrax DMR uses <Group name="BEAM0000"> … </Group>
    for g in root.iter():
        # tag may be '{namespace}Group'
        if g.tag.endswith('Group') and 'name' in g.attrib:
            name = g.attrib['name']
            if name.upper().startswith("BEAM"):
                beams.add(name)
    return sorted(beams)

def discover_beam_var_paths(dmr_xml_text: str) -> dict:
    """
    Parse DMR XML and return:
      { 'BEAM0000': {'agbd':'agbd', 'lat_lowestmode':'geolocation/lat_lowestmode', ...}, ... }
    Paths are relative to the BEAM group and include subgroup prefixes when present.
    """
    root = ET.fromstring(dmr_xml_text)

    def local(tag):  # strip namespace
        return tag.split('}')[-1]

    def is_group(e): return local(e.tag) == 'Group'

    def is_variable_elem(e):
        # In DAP4 DMR, variables are typed elements (Float32, Int64, String, Structure, etc.)
        # Anything with a 'name' that is NOT a Group/Attribute/Dimension/EnumDef/Map is a variable.
        l = local(e.tag)
        if l in ('Group', 'Attribute', 'Dimension', 'EnumDef', 'Map'):
            return False
        return 'name' in e.attrib

    beam_maps = {}

    def walk(group_elem, prefix, mapping):
        # record variables at this level
        for child in list(group_elem):
            if is_variable_elem(child):
                vname = child.attrib.get('name', '')
                if vname:
                    # keep first (shortest) path we see
                    mapping.setdefault(vname, f"{prefix}{vname}" if prefix else vname)
        # recurse into child groups
        for child in list(group_elem):
            if is_group(child):
                gname = child.attrib.get('name', '')
                if gname:
                    walk(child, f"{prefix}{gname}/" if prefix else f"{gname}/", mapping)

    # find each BEAM group and collect its variables (with relative paths)
    for elem in root.iter():
        if is_group(elem) and elem.attrib.get('name', '').upper().startswith('BEAM'):
            beam = elem.attrib['name']
            m = {}
            walk(elem, '', m)
            beam_maps[beam] = m

    return beam_maps


def build_ce_from_var_paths(beam_maps: dict):
    """
    Return (list_of_paths, beams_used).
    Each path is like '/BEAM0000/agbd' or '/BEAM0000/geolocation/lat_lowestmode'.
    """
    want = {
        'agbd': ['agbd'],
        'agbd_se': ['agbd_se'],
        'l4_quality_flag': ['l4_quality_flag'],
        'degrade_flag': ['degrade_flag'],
        'shot_number': ['shot_number'],
        'lat': ['lat_lowestmode', 'latitude'],
        'lon': ['lon_lowestmode', 'longitude'],
    }

    paths, beams_used = [], []

    for beam, name_to_path in beam_maps.items():
        lat_path = next((name_to_path[n] for n in want['lat'] if n in name_to_path), None)
        lon_path = next((name_to_path[n] for n in want['lon'] if n in name_to_path), None)
        if 'agbd' not in name_to_path or lat_path is None or lon_path is None:
            continue

        beams_used.append(beam)

        # core vars if present
        for k in ['agbd','agbd_se','l4_quality_flag','degrade_flag','shot_number']:
            if k in name_to_path:
                paths.append(f"/{beam}/{name_to_path[k]}")
        # required lat/lon
        paths.append(f"/{beam}/{lat_path}")
        paths.append(f"/{beam}/{lon_path}")

    return paths, beams_used

def build_dap4_ce(beams: List[str]) -> str:
    """
    Build a DAP4 constraint listing only variables we need per beam.
    """
    want = ["agbd","agbd_se","l4_quality_flag","degrade_flag","shot_number",
            "lat_lowestmode","lon_lowestmode","latitude","longitude"]
    parts=[]
    for b in beams:
        for v in want:
            parts.append(f"/{b}/{v}")
    return "?dap4.ce=" + ";".join(parts)

def read_gedi_from_nc_bytes(nc_bytes: bytes, beams: list, aoi_bbox):
    xmin, ymin, xmax, ymax = aoi_bbox
    rows=[]
    ds = None
    try:
        # try in-memory
        ds = nc.Dataset("inmemory", mode="r", memory=nc_bytes)
    except Exception:
        # fallback: write to temp file and open
        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
            tmp.write(nc_bytes)
            tmp_path = tmp.name
        try:
            ds = nc.Dataset(tmp_path, mode="r")
        finally:
            try:
                os.remove(tmp_path)
            except Exception:
                pass
    def get_arr(grp, *names):
        for n in names:
            if n in grp.variables: return grp.variables[n][:]
        # fallback to nested geolocation group if present
        if "geolocation" in grp.groups:
            g = grp.groups["geolocation"]
            for n in names:
                if n in g.variables: return g.variables[n][:]
        return None
    for b in beams:
        if b not in ds.groups:  # if Hyrax dropped empty beams in subset
            continue
        grp = ds.groups[b]
        agbd    = get_arr(grp, "agbd")
        agbd_se = get_arr(grp, "agbd_se")
        qflag   = get_arr(grp, "l4_quality_flag")
        dflag   = get_arr(grp, "degrade_flag")
        shot    = get_arr(grp, "shot_number")
        lat     = get_arr(grp, "lat_lowestmode", "latitude")
        lon     = get_arr(grp, "lon_lowestmode", "longitude")
        if any(v is None for v in [agbd, agbd_se, qflag, dflag, shot, lat, lon]):
            continue
        n = len(agbd)
        # align lengths if needed
        for arr in [agbd_se, qflag, dflag, shot, lat, lon]:
            if len(arr) != n:
                m = min([n, len(agbd_se), len(qflag), len(dflag), len(shot), len(lat), len(lon)])
                agbd, agbd_se, qflag, dflag, shot, lat, lon = \
                    agbd[:m], agbd_se[:m], qflag[:m], dflag[:m], shot[:m], lat[:m], lon[:m]
                n = m
                break
        # filters
        good = (qflag == 1) & (dflag == 0)
        in_box = (lon >= xmin) & (lon <= xmax) & (lat >= ymin) & (lat <= ymax)
        mask = good & in_box & np.isfinite(agbd) & np.isfinite(lat) & np.isfinite(lon)
        if not np.any(mask): continue
        idx = np.where(mask)[0]
        for i in idx:
            rows.append({
                "footprint_id": f"{int(shot[i])}_{b}",
                "shot_number": int(shot[i]),
                "beam": b,
                "lat": float(lat[i]),
                "lon": float(lon[i]),
                "agbd": float(agbd[i]),
                "agbd_se": float(agbd_se[i]),
            })
    return pd.DataFrame(rows)

def read_gedi_from_h5(local_path: str, aoi_bbox, want_quality=True) -> pd.DataFrame:
    xmin, ymin, xmax, ymax = aoi_bbox
    rows = []
    with h5py.File(local_path, "r") as f:
        beams = [k for k in f.keys() if k.upper().startswith("BEAM")]
        for b in beams:
            grp = f[b]

            def get_arr(names):
                # try at beam root, then geolocation subgroup
                for n in names:
                    if n in grp:
                        return grp[n][()]
                if "geolocation" in grp:
                    g = grp["geolocation"]
                    for n in names:
                        if n in g:
                            return g[n][()]
                return None

            agbd    = get_arr(["agbd"])
            agbd_se = get_arr(["agbd_se"])
            qflag   = get_arr(["l4_quality_flag"])
            dflag   = get_arr(["degrade_flag"])
            shot    = get_arr(["shot_number"])
            lat     = get_arr(["lat_lowestmode", "latitude"])
            lon     = get_arr(["lon_lowestmode", "longitude"])

            if any(v is None for v in [agbd, agbd_se, qflag, dflag, shot, lat, lon]):
                continue

            # ensure 1D and aligned
            n = min(len(agbd), len(agbd_se), len(qflag), len(dflag), len(shot), len(lat), len(lon))
            agbd, agbd_se, qflag, dflag, shot, lat, lon = \
                agbd[:n], agbd_se[:n], qflag[:n], dflag[:n], shot[:n], lat[:n], lon[:n]

            mask = np.isfinite(agbd) & np.isfinite(lat) & np.isfinite(lon)
            if want_quality:
                mask &= (qflag == 1) & (dflag == 0)
            mask &= (lon >= xmin) & (lon <= xmax) & (lat >= ymin) & (lat <= ymax)
            idx = np.where(mask)[0]
            if idx.size == 0:
                continue

            for i in idx:
                rows.append({
                    "footprint_id": f"{int(shot[i])}_{b}",
                    "shot_number": int(shot[i]),
                    "beam": b,
                    "lat": float(lat[i]),
                    "lon": float(lon[i]),
                    "agbd": float(agbd[i]),
                    "agbd_se": float(agbd_se[i]),
                })
    return pd.DataFrame(rows)

# ---------------- HLS helpers ----------------
def stac_search_hls(aoi_bbox, start_dt, end_dt, max_items=MAX_HLS_ITEMS_PER_MONTH):
    client = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
    search = client.search(
        collections=["HLSL30.v2.0","HLSS30.v2.0"],
        bbox=list(aoi_bbox),
        datetime=f"{start_dt.strftime('%Y-%m-%d')}/{end_dt.strftime('%Y-%m-%d')}",
        limit=10000
    )
    items = list(search.get_items())
    if not items: return []
    items.sort(key=lambda it: it.datetime or dtparser.parse(it.properties.get("datetime")))
    if len(items) > max_items: items = items[:max_items]
    return [pc.sign(it) for it in items]

def item_band_map(item):
    coll = (item.collection_id or item.to_dict().get("collection", ""))
    if "HLSS30" in coll:
        return {"B02":"B02","B03":"B03","B04":"B04","B05":"B08","B06":"B11","B07":"B12","Fmask":"Fmask"}
    else:
        return {"B02":"B02","B03":"B03","B04":"B04","B05":"B05","B06":"B06","B07":"B07","Fmask":"Fmask"}

def sample_items_at_points(items, pts_lonlat, allowed_fmask=FMASK_CLEAR_VALUES):
    n = len(pts_lonlat)
    per_band_values = {b: [] for b in ["B02","B03","B04","B05","B06","B07"]}
    gdal_env = rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="10000000"
    )
    with gdal_env:
        for it in items:
            bmap = item_band_map(it)
            fmask_asset = it.assets.get(bmap["Fmask"])
            if fmask_asset is None: continue
            try:
                with rasterio.open(fmask_asset.href) as src_m:
                    fmask_vals = np.array([v[0] for v in src_m.sample(pts_lonlat)], dtype=np.float32)
            except RasterioIOError:
                continue
            keep = np.isin(fmask_vals, list(allowed_fmask))
            if not keep.any(): continue
            for band in ["B02","B03","B04","B05","B06","B07"]:
                asset = it.assets.get(bmap[band])
                if asset is None:
                    per_band_values[band].append(np.full(n, np.nan, dtype=np.float32)); continue
                try:
                    with rasterio.open(asset.href) as src_b:
                        vals = np.array([v[0] for v in src_b.sample(pts_lonlat)], dtype=np.float32)
                except RasterioIOError:
                    vals = np.full(n, np.nan, dtype=np.float32)
                vals[~keep] = np.nan
                per_band_values[band].append(vals)
    return per_band_values

def nanmedian_stack(arr_list):
    if not arr_list: return None
    stack = np.vstack(arr_list)
    return np.nanmedian(stack, axis=0)


# --- QA decoding for HLS v2.0 Fmask (bit-packed) ---
# Bits (MSB→LSB): 7-6 aerosol, 5 water, 4 snow/ice, 3 shadow, 2 adjacent, 1 cloud, 0 (reserved)
# Ref: HLS User Guide v2.0, Table 9.  (HLSS30/ HLSL30 v2.0)  # ← documentation note

# Configure what to allow:
HLS_QA_ACCEPT_WATER    = False   # set True if you want to keep clear water pixels
HLS_QA_ACCEPT_ADJACENT = False   # set True to keep pixels adjacent to cloud/shadow
HLS_QA_ACCEPT_SNOW     = False   # set True to keep snow/ice

def hls_v2_clear_mask(qa_array: np.ndarray,
                      accept_water=False,
                      accept_adjacent=False,
                      accept_snow=False) -> np.ndarray:
    qa = qa_array.astype(np.uint16)
    BAD = 0
    BAD |= (1 << 1)  # cloud
    BAD |= (1 << 3)  # cloud shadow
    if not accept_adjacent:
        BAD |= (1 << 2)  # adjacent to cloud/shadow
    if not accept_snow:
        BAD |= (1 << 4)  # snow/ice
    if not accept_water:
        BAD |= (1 << 5)  # water
    # Ignore aerosol level bits (6–7) and reserved bit 0
    return (qa & BAD) == 0


# ---------------- Main ----------------
def main(args):
    # --- Auth ---
    earthaccess.login(strategy="netrc")

    # --- Resolve config from CLI ---
    AOI_BBOX = _resolve_aoi_bbox(args)
    YEAR = args.year
    MONTHS = _parse_months(args.months)
    ROLLING_DAYS = args.rolling_days
    MAX_HLS_ITEMS_PER_MONTH = args.max_hls_per_month
    GEDI_L4A_COLLECTION_ID = args.gedi_collection_id

    # HLS QA accept flags (used by hls_v2_clear_mask)
    global HLS_QA_ACCEPT_WATER, HLS_QA_ACCEPT_ADJACENT, HLS_QA_ACCEPT_SNOW
    HLS_QA_ACCEPT_WATER = bool(args.accept_water)
    HLS_QA_ACCEPT_ADJACENT = bool(args.accept_adjacent)
    HLS_QA_ACCEPT_SNOW = bool(args.accept_snow)

    # Workspace & outputs
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)
    gedi_cache = os.path.join(workdir, "gedi_cache")
    hls_cache  = os.path.join(workdir, "hls_cache_hls_v2")
    os.makedirs(gedi_cache, exist_ok=True)
    os.makedirs(hls_cache,  exist_ok=True)

    OUT_CSV = args.out_csv if os.path.isabs(args.out_csv) else os.path.join(workdir, args.out_csv)

    # GEDI per-month cap
    max_per_month = args.limit_gedi_per_month if args.limit_gedi_per_month is not None else None

    print(f"AOI bbox: {AOI_BBOX}")
    print(f"Year: {YEAR}  Months: {MONTHS}  Rolling ±{ROLLING_DAYS}d")
    print(f"Workdir: {workdir}")
    print(f"GEDI cache: {gedi_cache} | HLS cache: {hls_cache}")
    print(f"Output CSV: {OUT_CSV}")

    # -----------------------------
    # 1) GEDI L4A: download & read
    # -----------------------------
    gedi_rows = []
    print(f"Downloading + reading GEDI granules for months {sorted(set(MONTHS))}…")

    for m in sorted(set(MONTHS)):
        m_start = datetime(YEAR, m, 1)
        m_end   = (m_start + relativedelta(months=1)) - timedelta(days=1)
        q_start = m_start - timedelta(days=ROLLING_DAYS)
        q_end   = m_end + timedelta(days=ROLLING_DAYS)

        # Search via concept_id; fallback to short_name
        try:
            gedi_results = earthaccess.search_data(
                concept_id=GEDI_L4A_COLLECTION_ID,
                bounding_box=(AOI_BBOX[0], AOI_BBOX[1], AOI_BBOX[2], AOI_BBOX[3]),
                temporal=(q_start, q_end)
            )
        except Exception:
            gedi_results = earthaccess.search_data(
                short_name="GEDI_L4A_AGB_Density_V2_1",
                bounding_box=(AOI_BBOX[0], AOI_BBOX[1], AOI_BBOX[2], AOI_BBOX[3]),
                temporal=(q_start, q_end)
            )

        if not gedi_results:
            print(f"  {YEAR}-{m:02d}: 0 GEDI granules")
            continue

        if max_per_month:
            gedi_results = gedi_results[:max_per_month]

        print(f"  {YEAR}-{m:02d}: {len(gedi_results)} GEDI granules")

        # Batch download (earthaccess skips existing files automatically)
        try:
            gedi_files = earthaccess.download(gedi_results, local_path=gedi_cache)
        except TypeError:
            gedi_files = earthaccess.download(gedi_results, gedi_cache)

        local_h5s = [str(p) for p in gedi_files if str(p).lower().endswith(".h5")]
        for fp in tqdm(local_h5s, desc=f"GEDI {YEAR}-{m:02d}", leave=False):
            try:
                df = read_gedi_from_h5(fp, AOI_BBOX, want_quality=True)
                if not df.empty:
                    gedi_rows.append(df)
            except Exception as e:
                warnings.warn(f"GEDI read failed for {fp}: {e}")

    if not gedi_rows:
        sys.exit("No GEDI shots passed AOI + QC filters. Try expanding AOI or time window.")

    gedi = pd.concat(gedi_rows, ignore_index=True).drop_duplicates(subset=["footprint_id"])
    print(f"GEDI shots kept: {len(gedi)}")

    # --------------------------------------------
    # 2) HLSS30 v2.0 (Earthdata): sample + median
    # --------------------------------------------
    # Nested helpers so you don't have to edit imports elsewhere.
    from rasterio.warp import transform as rio_transform
    from rasterio.crs import CRS

    def search_hlss30(aoi_bbox, start_dt, end_dt, limit=None):
        """Search Earthdata CMR for BOTH HLSS30 and HLSL30 v2.0 granules."""

        def _q(short_name):
            try:
                return earthaccess.search_data(
                    short_name=short_name,
                    version="2.0",
                    bounding_box=(aoi_bbox[0], aoi_bbox[1], aoi_bbox[2], aoi_bbox[3]),
                    temporal=(start_dt, end_dt),
                )
            except Exception:
                return []

        s30 = _q("HLSS30")
        l30 = _q("HLSL30")
        results = list(s30) + list(l30)

        # sort by start time if available
        try:
            results = sorted(
                results,
                key=lambda g: pd.to_datetime(
                    getattr(g, "umm", {})
                    .get("TemporalExtent", {})
                    .get("RangeDateTime", {})
                    .get("BeginningDateTime", "1970-01-01")
                ),
            )
        except Exception:
            pass

        if limit:
            results = results[:limit]
        return results

    def band_paths_from_downloaded(files):
        """
        Return canonical band paths for one HLS granule (HLSS30 or HLSL30).
        Maps:
          Blue=B02, Green=B03, Red=B04,
          NIR=B08 (S2) or B05 (L8/9),
          SWIR1=B11 (S2) or B06 (L8/9),
          SWIR2=B12 (S2) or B07 (L8/9)
        """
        paths = {
            b: None for b in ["Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2", "Fmask"]
        }
        for p in files:
            s = str(p)
            low = s.lower()
            if not low.endswith(".tif"):
                continue

            if re.search(r"\.b02\.tif$", low):
                paths["Blue"] = s
            elif re.search(r"\.b03\.tif$", low):
                paths["Green"] = s
            elif re.search(r"\.b04\.tif$", low):
                paths["Red"] = s

            # NIR
            elif re.search(r"\.b08\.tif$", low):
                paths["NIR"] = s  # Sentinel-2
            elif re.search(r"\.b05\.tif$", low) and paths["NIR"] is None:  # Landsat
                paths["NIR"] = s

            # SWIR1
            elif re.search(r"\.b11\.tif$", low):
                paths["SWIR1"] = s  # Sentinel-2
            elif re.search(r"\.b06\.tif$", low) and paths["SWIR1"] is None:  # Landsat
                paths["SWIR1"] = s

            # SWIR2
            elif re.search(r"\.b12\.tif$", low):
                paths["SWIR2"] = s  # Sentinel-2
            elif re.search(r"\.b07\.tif$", low) and paths["SWIR2"] is None:  # Landsat
                paths["SWIR2"] = s

            elif re.search(r"\.fmask\.tif$", low):
                paths["Fmask"] = s
        return paths

    def sample_hlss30_month(aoi_bbox, start_dt, end_dt, pts_lonlat, max_items):
        """
        Download HLS (HLSS30 + HLSL30) v2.0 granules from Earthdata, sample only the GEDI points
        that fall inside each tile (in tile CRS), mask by QA/Fmask, and return rolling-window medians.
        """
        # 1) Search BOTH sensors
        results = search_hlss30(aoi_bbox, start_dt, end_dt, limit=max_items)
        if not results:
            return None, 0

        # 2) Download assets to cache
        try:
            dl = earthaccess.download(results, local_path=hls_cache)
        except TypeError:
            dl = earthaccess.download(results, hls_cache)

        # 3) Bundle files by granule (strip trailing band suffix for either sensor)
        bundles = {}
        for p in dl:
            s = str(p)
            low = s.lower()
            if not low.endswith(".tif"):
                continue
            base = re.sub(
                r"\.(b02|b03|b04|b05|b06|b07|b08|b11|b12|fmask)\.tif$",
                "",
                s,
                flags=re.I,
            )
            bundles.setdefault(base, []).append(s)

        # 4) Prepare containers
        n = len(pts_lonlat)
        lons = np.array([xy[0] for xy in pts_lonlat], dtype=np.float64)
        lats = np.array([xy[1] for xy in pts_lonlat], dtype=np.float64)
        bands = ["Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"]
        per_band = {b: [] for b in bands}

        used = 0

        # 5) Process each granule
        for base, files in tqdm(bundles.items(), desc="HLS sampling", leave=False):
            paths = band_paths_from_downloaded(files)
            if paths["Fmask"] is None:
                continue  # need Fmask to decide 'clear' pixels

            try:
                with rasterio.open(paths["Fmask"]) as src_m:
                    xs, ys = rio_transform(CRS.from_epsg(4326), src_m.crs, lons, lats)
                    xs = np.asarray(xs)
                    ys = np.asarray(ys)
                    left, bottom, right, top = src_m.bounds
                    in_tile = (
                        (xs >= left) & (xs <= right) & (ys >= bottom) & (ys <= top)
                    )
                    idx_tile = np.where(in_tile)[0]
                    if idx_tile.size == 0:
                        continue

                    coords_tile = list(zip(xs[idx_tile], ys[idx_tile]))
                    fmask_vals = np.array(
                        [v[0] for v in src_m.sample(coords_tile)], dtype=np.uint16
                    )

                keep_local = hls_v2_clear_mask(
                    fmask_vals,
                    accept_water=HLS_QA_ACCEPT_WATER,
                    accept_adjacent=HLS_QA_ACCEPT_ADJACENT,
                    accept_snow=HLS_QA_ACCEPT_SNOW,
                )
                if not keep_local.any():
                    continue

                idx_keep = idx_tile[keep_local]

                # sample canonical bands
                for b in bands:
                    pth = paths[b]
                    full = np.full(n, np.nan, dtype=np.float32)
                    if pth is not None and idx_keep.size:
                        with rasterio.open(pth) as src_b:
                            vals = np.array(
                                [
                                    v[0]
                                    for v in src_b.sample(
                                        list(zip(xs[idx_keep], ys[idx_keep]))
                                    )
                                ],
                                dtype=np.float32,
                            )
                        full[idx_keep] = vals
                    per_band[b].append(full)

                used += 1

            except Exception as e:
                warnings.warn(f"HLS sampling failed for {base}: {e}")
                continue

        if used == 0:
            return None, 0

        # 6) Median across the stack (both sensors together)
        med = {}
        for k, stacks in per_band.items():
            if len(stacks) == 0:
                med[k] = np.full(n, np.nan, dtype=np.float32)
            else:
                med[k] = np.nanmedian(np.vstack(stacks), axis=0)

        return med, used

    # Prepare GEDI points once
    pts = list(zip(gedi["lon"].to_numpy(), gedi["lat"].to_numpy()))
    all_months = []

    for m in MONTHS:
        m_start = datetime(YEAR, m, 1)
        m_end   = (m_start + relativedelta(months=1)) - timedelta(days=1)
        q_start = m_start - timedelta(days=ROLLING_DAYS)
        q_end   = m_end + timedelta(days=ROLLING_DAYS)

        med, used = sample_hlss30_month(AOI_BBOX, q_start, q_end, pts, MAX_HLS_ITEMS_PER_MONTH)
        if med is None:
            print(f"[warn] No HLSS30 items for {YEAR}-{m:02d}")
            out = gedi[["footprint_id","lon","lat","agbd","agbd_se","beam","shot_number"]].copy()
            out["ym"] = f"{YEAR}-{m:02d}"
            for b in ["Blue","Green","Red","NIR","SWIR1","SWIR2"]:
                out[b] = np.nan
            all_months.append(out)
            continue

        out = gedi[["footprint_id","lon","lat","agbd","agbd_se","beam","shot_number"]].copy()
        out["ym"]   = f"{YEAR}-{m:02d}"
        out["Blue"] = med["Blue"]
        out["Green"] = med["Green"]
        out["Red"] = med["Red"]
        out["NIR"] = med["NIR"]
        out["SWIR1"] = med["SWIR1"]
        out["SWIR2"] = med["SWIR2"]
        print(f"[ok] {YEAR}-{m:02d}: HLS items used={used}, rows={len(out)}")

        all_months.append(out)

    # --------------------------------
    # 3) Save combined CSV
    # --------------------------------
    final = pd.concat(all_months, ignore_index=True)
    final.to_csv(OUT_CSV, index=False)
    print(f"Saved: {OUT_CSV}")
    print(final.head())



if __name__ == "__main__":
    args = parse_cli()
    main(args)
