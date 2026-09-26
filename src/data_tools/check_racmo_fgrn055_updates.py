#!/usr/bin/env python3
"""
check_racmo_fgrn055_updates.py

Checks for updates to the RACMO2.3p2 FGRN055 (Greenland, 5.5 km) surface
mass balance / climate dataset, using the two sources the group actually
watches:

  1. https://www.projects.science.uu.nl/iceclimate/models/racmo-data.php#4-1
     — IMAU's own "Data Greenland" table, which lists every RACMO
     model/domain combination for Greenland (currently RACMO2.3p2 FGRN055,
     RACMO2.4p1 FGRN055, RACMO2.4p1 FGRN11, RACMO2.3p3 FGRN11, plus the
     historic/21st-century CESM2-forced FGRN11 runs) along with each one's
     resolution, period, forcing, and availability status.
  2. https://zenodo.org/records/19254871 ("Monthly averaged RACMO2.3p2
     variables: Greenland") via the Zenodo REST API — the actual archive
     of monthly NetCDF files for RACMO2.3p2 FGRN055.

This script reports two things the group cares about:

  A. A VERSION UPGRADE: a Greenland RACMO product other than the current
     RACMO2.3p2 FGRN055 (e.g. RACMO2.4p1 FGRN055/FGRN11, RACMO2.3p3 FGRN11,
     or a "1 km downscaled" product) that is at 5.5 km resolution or finer
     AND is actually available (not "in preparation" / "not operational").
     This is tracked by diffing IMAU's Greenland table between runs, since
     that page is where such a switch (e.g. RACMO2.4p1 FGRN055 flipping
     from "Available, but not operational" to a live Zenodo link) shows up
     first.
  B. NEW MONTHS/YEARS in the existing RACMO2.3p2 FGRN055 Zenodo archive:
     each variable is one NetCDF file per record, named with an embedded
     date range (e.g. "..._193909_202512.nc"); an advancing end-date means
     new months were appended. A new Zenodo version number, or a changed
     file checksum/size for an unchanged date range, means the existing
     period was reprocessed/corrected rather than merely extended.

State is cached locally in one JSON file (default: racmo_fgrn055_state.json
next to this script; override with --state-file). By default this script
is READ-ONLY: it always fetches the current state and, if a state file
already exists, diffs and reports against it — but it never writes or
overwrites that file unless you pass --update.

Usage:
    python3 check_racmo_fgrn055_updates.py                          # report only, no file written
    python3 check_racmo_fgrn055_updates.py --update                 # report, then save/update the state file
    python3 check_racmo_fgrn055_updates.py --min-resolution-km 5.5   # threshold for "meets your criterion" (default)
    python3 check_racmo_fgrn055_updates.py --json                    # also print the fetched state as raw JSON

No third-party dependencies required (uses urllib from the standard library).
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

UU_PAGE_URL = "https://www.projects.science.uu.nl/iceclimate/models/racmo-data.php"
ZENODO_API_BASE = "https://zenodo.org/api"
DEFAULT_ZENODO_RECORD = "19254871"  # "Monthly averaged RACMO2.3p2 variables: Greenland", v1
DEFAULT_STATE_FILE = Path(__file__).resolve().parent / "racmo_fgrn055_state.json"
DEFAULT_MIN_RESOLUTION_KM = 5.5
TARGET_CATEGORY = "Data Greenland"
TIMEOUT = 60
USER_AGENT = "racmo-fgrn055-update-checker/1.0"


# --------------------------------------------------------------------------
# HTTP helpers
# --------------------------------------------------------------------------

def _http_get_json(url: str, params: dict | None = None) -> dict:
    full_url = f"{url}?{urllib.parse.urlencode(params)}" if params else url
    req = urllib.request.Request(full_url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Request failed ({e.code}) for {full_url}: {e.read().decode(errors='replace')}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach {full_url}: {e}") from e


def _http_get_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Request failed ({e.code}) for {url}: {e.read().decode(errors='replace')}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach {url}: {e}") from e


# --------------------------------------------------------------------------
# IMAU "Data Greenland" table
# --------------------------------------------------------------------------

_TOKEN_RE = re.compile(
    r'<h2[^>]*>(?P<cat>Data [^<]+)</h2>'
    r'|<h5>(?P<sec>[^<]+)</h5>'
    r'|<table style="margin-top:10px">(?P<entry>.*?)</table>',
    re.DOTALL,
)
_FIELD_RE = re.compile(r'<h4>\s*([^<:]+):\s*</h4>.*?<p>(.*?)</p>', re.DOTALL)
_LINK_RE = re.compile(r'href="([^"]+)"')
_RES_KM_RE = re.compile(r'(\d+(?:\.\d+)?)\.?\s*km', re.IGNORECASE)
_ONE_KM_RE = re.compile(r'\b1\s*km\b', re.IGNORECASE)


def _clean(s: str) -> str:
    s = re.sub(r'<[^>]+>', ' ', s)
    s = html.unescape(s)
    s = re.sub(r'\s+', ' ', s).strip()
    s = re.sub(r'\s+([,.;])', r'\1', s)
    return s


def fetch_uu_entries(page_url: str, category_filter: str | None) -> dict:
    """Parse every dataset "card" on the IMAU RACMO data page into a dict keyed by
    "<category> | <section> | <name>", optionally restricted to one category
    (e.g. "Data Greenland").
    """
    text = _http_get_text(page_url)

    cat, sec = None, None
    entries: dict[str, dict] = {}
    for m in _TOKEN_RE.finditer(text):
        if m.group('cat'):
            cat = _clean(m.group('cat'))
        elif m.group('sec'):
            sec = _clean(m.group('sec'))
        elif m.group('entry') is not None:
            if category_filter and cat != category_filter:
                continue
            block = m.group('entry')
            fields = {_clean(k).lower(): _clean(v) for k, v in _FIELD_RE.findall(block)}
            links = _LINK_RE.findall(block)
            name = fields.get('name', '(unnamed)')
            key = f"{cat} | {sec} | {name}"

            resolution_text = fields.get('resolution', '')
            availability_text = fields.get('racmo data', '')
            res_match = _RES_KM_RE.search(resolution_text)
            resolution_km = float(res_match.group(1)) if res_match else None
            has_1km_downscaled = bool(_ONE_KM_RE.search(resolution_text + ' ' + availability_text))
            is_operational = not re.search(r'not operational|in preparation', availability_text, re.IGNORECASE)

            entries[key] = {
                "category": cat,
                "section": sec,
                "name": name,
                "resolution_text": resolution_text,
                "resolution_km": resolution_km,
                "has_1km_downscaled": has_1km_downscaled,
                "domain": fields.get('domain', ''),
                "period": fields.get('period', ''),
                "forcing": fields.get('forcing', ''),
                "availability_text": availability_text,
                "is_operational": is_operational,
                "links": links,
            }
    return entries


# --------------------------------------------------------------------------
# Zenodo record
# --------------------------------------------------------------------------

_DATE_RANGE_RE = re.compile(r'_(\d{6})_(\d{6})\.nc$')


def _latest_zenodo_record(requested_record_id: str) -> dict:
    """Given any record id in a Zenodo version series, return the JSON for
    whichever record is currently the latest version of that series.
    """
    base_record = _http_get_json(f"{ZENODO_API_BASE}/records/{requested_record_id}")
    concept_id = base_record.get("conceptrecid")
    if not concept_id:
        return base_record

    search = _http_get_json(
        f"{ZENODO_API_BASE}/records",
        {"q": f"conceptrecid:{concept_id}", "all_versions": "true", "sort": "-version", "size": 5},
    )
    hits = search.get("hits", {}).get("hits", [])
    if not hits:
        return base_record
    return hits[0]


def fetch_zenodo_state(requested_record_id: str) -> dict:
    record = _latest_zenodo_record(requested_record_id)
    md = record.get("metadata", {})

    files = {}
    coverage_start, coverage_end = None, None
    for f in record.get("files", []):
        key = f.get("key", "unknown")
        files[key] = {"size": f.get("size"), "checksum": f.get("checksum")}
        m = _DATE_RANGE_RE.search(key)
        if m:
            start, end = m.group(1), m.group(2)
            coverage_start = start if coverage_start is None or start < coverage_start else coverage_start
            coverage_end = end if coverage_end is None or end > coverage_end else coverage_end

    return {
        "requested_record_id": requested_record_id,
        "conceptrecid": record.get("conceptrecid"),
        "conceptdoi": record.get("conceptdoi"),
        "latest_record_id": record.get("id"),
        "latest_doi": record.get("doi"),
        "latest_version_label": md.get("version"),
        "title": md.get("title"),
        "publication_date": md.get("publication_date"),
        "updated": record.get("updated"),
        "coverage_start_yyyymm": coverage_start,
        "coverage_end_yyyymm": coverage_end,
        "files": files,
    }


# --------------------------------------------------------------------------
# Combined state
# --------------------------------------------------------------------------

def build_state(zenodo_record_id: str, page_url: str) -> dict:
    return {
        "checked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "zenodo": fetch_zenodo_state(zenodo_record_id),
        "uu_page": {
            "source_url": page_url,
            "entries": fetch_uu_entries(page_url, TARGET_CATEGORY),
        },
    }


# --------------------------------------------------------------------------
# Diffing and reporting
# --------------------------------------------------------------------------

def diff_zenodo(old: dict, new: dict) -> dict:
    old_files, new_files = old["files"], new["files"]
    return {
        "version_changed": old.get("latest_record_id") != new.get("latest_record_id"),
        "coverage_start_changed": old.get("coverage_start_yyyymm") != new.get("coverage_start_yyyymm"),
        "coverage_end_changed": old.get("coverage_end_yyyymm") != new.get("coverage_end_yyyymm"),
        "added_files": sorted(set(new_files) - set(old_files)),
        "removed_files": sorted(set(old_files) - set(new_files)),
        "modified_files": sorted(
            k for k in (set(old_files) & set(new_files))
            if old_files[k].get("checksum") != new_files[k].get("checksum")
            or old_files[k].get("size") != new_files[k].get("size")
        ),
    }


def diff_uu_entries(old: dict, new: dict) -> dict:
    added = sorted(set(new) - set(old))
    removed = sorted(set(old) - set(new))
    changed = {}
    for key in sorted(set(old) & set(new)):
        o, n = old[key], new[key]
        deltas = {}
        for field in ("resolution_text", "period", "availability_text", "is_operational", "links"):
            if o.get(field) != n.get(field):
                deltas[field] = {"old": o.get(field), "new": n.get(field)}
        if deltas:
            changed[key] = deltas
    return {"added": added, "removed": removed, "changed": changed}


def notable_candidates(entries: dict, min_resolution_km: float) -> list[dict]:
    """Greenland RACMO products, other than the baseline RACMO2.3p2 FGRN055,
    that are both operational/available AND at min_resolution_km or finer
    (or explicitly noted as having a 1 km downscaled product) — i.e. the
    kind of entry that would answer "is there a higher-resolution successor
    to RACMO2.3p2 FGRN055 yet?"
    """
    out = []
    for key, e in entries.items():
        if e["name"] == "RACMO2.3p2 FGRN055":
            continue
        meets_res = (e["resolution_km"] is not None and e["resolution_km"] <= min_resolution_km) or e["has_1km_downscaled"]
        if meets_res and e["is_operational"]:
            out.append(e)
    return out


def format_report(old: dict | None, new: dict, min_resolution_km: float, update: bool) -> str:
    lines = []
    lines.append("RACMO2.3p2 FGRN055 (Greenland, 5.5 km) — dataset update check")
    lines.append("=" * 63)
    lines.append(f"Checked: {new['checked_at']}")
    lines.append("")

    z = new["zenodo"]
    lines.append(
        f"Zenodo archive: \"{z['title']}\" — record {z['latest_record_id']} "
        f"(version {z['latest_version_label'] or '1'}), DOI {z['latest_doi']}, "
        f"published {z['publication_date']}."
    )
    if z["coverage_start_yyyymm"] and z["coverage_end_yyyymm"]:
        lines.append(f"Data currently spans {z['coverage_start_yyyymm']} to {z['coverage_end_yyyymm']} (YYYYMM), across {len(z['files'])} file(s).")

    lines.append("")
    candidates = notable_candidates(new["uu_page"]["entries"], min_resolution_km)
    if candidates:
        lines.append(f"NOTABLE: Greenland RACMO product(s) at <= {min_resolution_km} km (or with a 1 km downscaled "
                      f"version) other than RACMO2.3p2 FGRN055 are currently available:")
        for c in candidates:
            lines.append(f"  - {c['name']}: {c['resolution_text']} | period {c['period']} | {c['availability_text']}")
    else:
        lines.append(
            f"No Greenland RACMO product at <= {min_resolution_km} km resolution other than RACMO2.3p2 FGRN055 "
            f"is currently marked available (checked RACMO2.4p1 FGRN055/FGRN11, RACMO2.3p3 FGRN11)."
        )

    lines.append("")
    if old is None:
        lines.append(
            "No prior check found" + (" — this run establishes the baseline." if update else ", and no state "
            "file was written (pass --update to save this as a baseline).") +
            " Re-run later to detect version changes, new months/years, or reprocessing."
        )
        return "\n".join(lines)

    lines.append(f"Comparing against previous check: {old['checked_at']}")
    lines.append("-" * 40)

    zd = diff_zenodo(old["zenodo"], new["zenodo"])
    ud = diff_uu_entries(old["uu_page"]["entries"], new["uu_page"]["entries"])

    anything_changed = (
        zd["version_changed"] or zd["coverage_start_changed"] or zd["coverage_end_changed"]
        or zd["added_files"] or zd["removed_files"] or zd["modified_files"]
        or ud["added"] or ud["removed"] or ud["changed"]
    )

    if not anything_changed:
        lines.append("No changes detected since the last check: same Zenodo version and file coverage, "
                      "same RACMO product listings on the IMAU Greenland data page.")
        return "\n".join(lines)

    if zd["version_changed"]:
        lines.append(
            f"* ZENODO VERSION CHANGE: record {old['zenodo']['latest_record_id']} "
            f"({old['zenodo']['publication_date']}) -> record {z['latest_record_id']} ({z['publication_date']})."
        )
    if zd["coverage_end_changed"]:
        lines.append(
            f"* NEW MONTHS ADDED: data end date advanced from {old['zenodo']['coverage_end_yyyymm']} "
            f"to {z['coverage_end_yyyymm']}."
        )
    if zd["coverage_start_changed"]:
        lines.append(
            f"* Data start date changed from {old['zenodo']['coverage_start_yyyymm']} to {z['coverage_start_yyyymm']}."
        )
    if zd["added_files"]:
        lines.append(f"* {len(zd['added_files'])} new file(s) in the Zenodo archive: {', '.join(zd['added_files'][:5])}"
                      + (" ..." if len(zd['added_files']) > 5 else ""))
    if zd["removed_files"]:
        lines.append(f"* {len(zd['removed_files'])} file(s) removed from the Zenodo archive: {', '.join(zd['removed_files'][:5])}"
                      + (" ..." if len(zd['removed_files']) > 5 else ""))
    if zd["modified_files"]:
        lines.append(f"* {len(zd['modified_files'])} file(s) reprocessed/replaced (checksum or size changed) "
                      f"without changing the covered date range: {', '.join(zd['modified_files'][:5])}"
                      + (" ..." if len(zd['modified_files']) > 5 else ""))

    for key in ud["added"]:
        e = new["uu_page"]["entries"][key]
        flag = " <-- meets your resolution criterion" if e in candidates else ""
        lines.append(f"* NEW ENTRY on IMAU Greenland data page: {e['name']} ({e['resolution_text']}, "
                      f"{e['availability_text']}){flag}")
    for key in ud["removed"]:
        lines.append(f"* ENTRY REMOVED from IMAU Greenland data page: {key.split(' | ')[-1]}")
    for key, deltas in ud["changed"].items():
        name = key.split(" | ")[-1]
        parts = [f"{field}: {d['old']!r} -> {d['new']!r}" for field, d in deltas.items()]
        lines.append(f"* CHANGED: {name} — " + "; ".join(parts))

    return "\n".join(lines)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def check_racmo(zenodo_record_id: str, page_url: str, state_file: Path, min_resolution_km: float,
                 update: bool = False) -> tuple[str, dict]:
    new_state = build_state(zenodo_record_id, page_url)

    old_state = None
    if state_file.exists():
        old_state = json.loads(state_file.read_text())

    report = format_report(old_state, new_state, min_resolution_km, update)

    if update:
        state_file.write_text(json.dumps(new_state, indent=2, sort_keys=True))
        report += f"\n\n(state saved to {state_file})"
    else:
        report += "\n\n(no state file written — pass --update to save this state)"

    return report, new_state


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--zenodo-record", default=DEFAULT_ZENODO_RECORD,
                         help=f"Zenodo record id for the FGRN055 archive (default: {DEFAULT_ZENODO_RECORD}); "
                              f"any version's id works, the latest version is resolved automatically")
    parser.add_argument("--page-url", default=UU_PAGE_URL,
                         help="IMAU RACMO data page URL to check for new/changed Greenland product listings")
    parser.add_argument("--min-resolution-km", type=float, default=DEFAULT_MIN_RESOLUTION_KM,
                         help=f"Flag Greenland RACMO products at or finer than this resolution (default: {DEFAULT_MIN_RESOLUTION_KM})")
    parser.add_argument("--state-file", type=Path, default=DEFAULT_STATE_FILE,
                         help="Where to cache state between runs (default: racmo_fgrn055_state.json next to this script)")
    parser.add_argument("--update", action="store_true",
                         help="Save the freshly fetched state to the state file (creating or overwriting it). "
                              "Default is read-only: report a diff (if a state file exists) but never write one.")
    parser.add_argument("--json", action="store_true", help="Also print the raw new-state JSON to stdout")
    args = parser.parse_args()

    try:
        report, new_state = check_racmo(
            args.zenodo_record, args.page_url, args.state_file, args.min_resolution_km, args.update
        )
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print(report)

    if args.json:
        print("\n--- raw state (JSON) ---")
        print(json.dumps(new_state, indent=2, sort_keys=True))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
