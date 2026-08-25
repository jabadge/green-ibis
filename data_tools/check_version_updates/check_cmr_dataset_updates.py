#!/usr/bin/env python3
"""
check_cmr_dataset_updates.py

Checks NASA's Common Metadata Repository (CMR) — the live catalog behind
every NASA-archived Earth science dataset, including all of NSIDC's for the 
current state of a dataset given its CMR "short name" (e.g.
NSIDC-0731 or IDBMG4) and reports:

  1. The latest released version of the dataset and when it was last updated.
  2. Whether any NEW years of data have been added since the last check.
  3. Whether any granules within YEARS/GRANULES ALREADY RELEASED have been
     replaced/reprocessed since the last check (detected via per-granule
     SHA-256 checksum, size, and CMR revision-id changes), even if the
     collection's headline version number did not change.

State is cached locally, one JSON file per dataset, named
"<short_name>_state.json" next to this script by default (override with
--state-file). By default this script is READ-ONLY: it always fetches the
current CMR state and, if a state file already exists, diffs and reports
against it — but it never writes or overwrites that file unless you pass
--update. Pass --update to save the freshly fetched state (creating the
file if it doesn't exist yet, or overwriting it after reporting the diff
if it does), so that a later run has something to compare against.

Usage:
    python3 check_cmr_dataset_updates.py NSIDC-0731        # report only, no file written
    python3 check_cmr_dataset_updates.py IDBMG4 --update   # report, then save/update the state file
    python3 check_cmr_dataset_updates.py NSIDC-0630 --state-file /path/to/state.json --update
    python3 check_cmr_dataset_updates.py NSIDC-0731 --json # also print the fetched state as raw JSON

No third-party dependencies required (uses urllib from the standard library).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

CMR_BASE = "https://cmr.earthdata.nasa.gov/search"
TIMEOUT = 60
USER_AGENT = "cmr-dataset-update-checker/1.0 (+https://cmr.earthdata.nasa.gov/search)"
SCRIPT_DIR = Path(__file__).resolve().parent


def default_state_file(short_name: str) -> Path:
    """One state file per dataset, named after its short_name, next to this script."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", short_name)
    return SCRIPT_DIR / f"{safe}_state.json"


# --------------------------------------------------------------------------
# CMR access helpers
# --------------------------------------------------------------------------

def _http_get(url: str, params: dict, extra_headers: dict | None = None):
    query = urllib.parse.urlencode(params)
    full_url = f"{url}?{query}"
    req = urllib.request.Request(full_url, headers={"User-Agent": USER_AGENT, **(extra_headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            body = json.loads(resp.read().decode("utf-8"))
            return body, resp.headers
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"CMR request failed ({e.code}) for {full_url}: {e.read().decode(errors='replace')}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach CMR at {full_url}: {e}") from e


def get_collections(short_name: str) -> list[dict]:
    """Every collection (dataset version) CMR currently has registered for short_name."""
    body, _ = _http_get(f"{CMR_BASE}/collections.umm_json", {"short_name": short_name, "page_size": 100})
    return body.get("items", [])


def pick_active_collection(collections: list[dict]) -> dict:
    """Pick the entry representing the current, actively distributed version.

    Prefers CollectionProgress == COMPLETE (the DAAC's flag for the current,
    non-deprecated release); falls back to the highest version number if
    that flag is absent on every entry.
    """
    def version_num(entry: dict) -> float:
        v = entry["umm"].get("Version", "0")
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0

    complete = [c for c in collections if c["umm"].get("CollectionProgress") == "COMPLETE"]
    pool = complete or collections
    return max(pool, key=version_num)


def get_all_granules(concept_id: str) -> list[dict]:
    """Page through every granule in a collection via CMR's granules.umm_json endpoint."""
    granules: list[dict] = []
    search_after = None
    while True:
        headers = {"CMR-Search-After": search_after} if search_after else None
        body, resp_headers = _http_get(
            f"{CMR_BASE}/granules.umm_json",
            {"collection_concept_id": concept_id, "page_size": 2000},
            headers,
        )
        page = body.get("items", [])
        granules.extend(page)
        search_after = resp_headers.get("CMR-Search-After")
        if not search_after or not page:
            break
    return granules


# --------------------------------------------------------------------------
# Summarizing CMR responses into a compact, diffable state
# --------------------------------------------------------------------------

def _granule_checksum(umm_granule: dict) -> str | None:
    for info in umm_granule.get("umm", {}).get("DataGranule", {}).get("ArchiveAndDistributionInformation", []):
        checksum = info.get("Checksum")
        if checksum:
            return f"{checksum.get('Algorithm', '?')}:{checksum.get('Value', '?')}"
    return None


def _granule_size_bytes(umm_granule: dict) -> int:
    total = 0
    for info in umm_granule.get("umm", {}).get("DataGranule", {}).get("ArchiveAndDistributionInformation", []):
        total += int(info.get("SizeInBytes", 0) or 0)
    return total


def _granule_year(umm_granule: dict) -> str | None:
    temporal = umm_granule.get("umm", {}).get("TemporalExtent", {}).get("RangeDateTime", {})
    begin = temporal.get("BeginningDateTime") or temporal.get("SingleDateTime")
    if not begin:
        single = umm_granule.get("umm", {}).get("TemporalExtent", {}).get("SingleDateTime")
        begin = single
    if not begin:
        return None
    return begin[:4]


def build_state(short_name: str) -> dict:
    collections = get_collections(short_name)
    if not collections:
        raise RuntimeError(f"No collections found in CMR for short_name={short_name!r}. "
                            f"Check the dataset ID (CMR 'short name') or CMR availability.")

    active = pick_active_collection(collections)
    active_umm = active["umm"]
    active_meta = active["meta"]

    granules = get_all_granules(active_meta["concept-id"])

    years: dict[str, dict] = {}
    for g in granules:
        year = _granule_year(g) or "unknown"
        gran_ur = g["umm"].get("GranuleUR", g["meta"].get("native-id", "unknown"))
        entry = {
            "revision_id": g["meta"].get("revision-id"),
            "revision_date": g["meta"].get("revision-date"),
            "production_date": g["umm"].get("DataGranule", {}).get("ProductionDateTime"),
            "checksum": _granule_checksum(g),
            "size_bytes": _granule_size_bytes(g),
        }
        year_bucket = years.setdefault(year, {"granules": {}})
        year_bucket["granules"][gran_ur] = entry

    for year_bucket in years.values():
        gset = year_bucket["granules"]
        year_bucket["granule_count"] = len(gset)
        year_bucket["total_size_bytes"] = sum(v["size_bytes"] for v in gset.values())
        latest = max((v["revision_date"] for v in gset.values() if v["revision_date"]), default=None)
        year_bucket["last_updated"] = latest

    doi = None
    doi_field = active_umm.get("DOI")
    if isinstance(doi_field, dict):
        doi = doi_field.get("DOI")

    temporal = active_umm.get("TemporalExtents", [{}])
    range_dt = temporal[0].get("RangeDateTimes", [{}])[0] if temporal else {}

    all_versions = []
    for c in collections:
        all_versions.append({
            "version": c["umm"].get("Version"),
            "concept_id": c["meta"]["concept-id"],
            "progress": c["umm"].get("CollectionProgress"),
            "revision_date": c["meta"].get("revision-date"),
        })
    all_versions.sort(key=lambda v: (v["version"] or ""))

    state = {
        "checked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "short_name": short_name,
        "entry_title": active_umm.get("EntryTitle"),
        "active_version": active_umm.get("Version"),
        "active_concept_id": active_meta["concept-id"],
        "active_collection_progress": active_umm.get("CollectionProgress"),
        "active_collection_updated": active_meta.get("revision-date"),
        "active_provider": active_meta.get("provider-id"),
        "doi": doi,
        "temporal_start": range_dt.get("BeginningDateTime"),
        "temporal_end": range_dt.get("EndingDateTime"),
        "landing_page": f"https://cmr.earthdata.nasa.gov/search/concepts/{active_meta['concept-id']}.html",
        "all_versions": all_versions,
        "years": years,
    }
    return state


# --------------------------------------------------------------------------
# Diffing two states and rendering a human-readable report
# --------------------------------------------------------------------------

def diff_states(old: dict, new: dict) -> dict:
    changes: dict = {
        "version_changed": old.get("active_version") != new.get("active_version"),
        "collection_metadata_changed": (
            old.get("active_version") == new.get("active_version")
            and old.get("active_collection_updated") != new.get("active_collection_updated")
        ),
        "new_years": sorted(set(new["years"]) - set(old["years"])),
        "removed_years": sorted(set(old["years"]) - set(new["years"])),
        "changed_years": {},  # year -> {added, removed, modified}
    }

    for year in sorted(set(old["years"]) & set(new["years"])):
        old_g = old["years"][year]["granules"]
        new_g = new["years"][year]["granules"]
        added = sorted(set(new_g) - set(old_g))
        removed = sorted(set(old_g) - set(new_g))
        modified = sorted(
            gid for gid in (set(old_g) & set(new_g))
            if old_g[gid].get("checksum") != new_g[gid].get("checksum")
            or old_g[gid].get("revision_id") != new_g[gid].get("revision_id")
            or old_g[gid].get("size_bytes") != new_g[gid].get("size_bytes")
        )
        if added or removed or modified:
            changes["changed_years"][year] = {"added": added, "removed": removed, "modified": modified}

    return changes


def format_report(old: dict | None, new: dict, changes: dict | None, update: bool) -> str:
    lines = []
    title = new.get("entry_title") or new["short_name"]
    header = f"{title} ({new['short_name']})"
    lines.append(header)
    lines.append("=" * len(header))
    lines.append(f"Checked: {new['checked_at']}")
    lines.append("")
    lines.append(
        f"Latest version: {new['active_version']} "
        f"(status: {new['active_collection_progress']}, provider: {new['active_provider']}), "
        f"last updated {new['active_collection_updated']}."
    )
    if new.get("doi"):
        lines.append(f"DOI: {new['doi']}")
    lines.append(f"CMR record: {new['landing_page']}")
    year_list = sorted(y for y in new["years"] if y != "unknown")
    coverage_note = f"{len(year_list)} year(s): {', '.join(year_list)}" if year_list else "no dated granules found"
    if "unknown" in new["years"]:
        coverage_note += " (plus granules with no parseable date)"
    lines.append(
        f"Temporal coverage: {new['temporal_start']} to {new['temporal_end']} ({coverage_note})."
    )

    lines.append("")
    if old is None or changes is None:
        if update:
            lines.append(
                "No prior check found — this run establishes the baseline. "
                "Re-run this script later (against the same --state-file) to detect "
                "new versions, new years of data, or revisions to already-released years."
            )
        else:
            lines.append(
                "No prior check found, and no state file was written (pass --update to save "
                "this as a baseline). Without a saved baseline, future runs have nothing to "
                "diff against."
            )
        return "\n".join(lines)

    lines.append(f"Comparing against previous check: {old['checked_at']}")
    lines.append("-" * 40)

    anything_changed = (
        changes["version_changed"]
        or changes["collection_metadata_changed"]
        or changes["new_years"]
        or changes["removed_years"]
        or changes["changed_years"]
    )

    if not anything_changed:
        lines.append("No changes detected since the last check: same version, same years, "
                      "same granule checksums/sizes/revisions throughout.")
        return "\n".join(lines)

    if changes["version_changed"]:
        lines.append(f"* VERSION CHANGE: {old['active_version']} -> {new['active_version']}.")

    if changes["collection_metadata_changed"]:
        lines.append(
            f"* Collection metadata for version {new['active_version']} was revised "
            f"({old['active_collection_updated']} -> {new['active_collection_updated']}) "
            f"without a version-number change."
        )

    if changes["new_years"]:
        for year in changes["new_years"]:
            n = new["years"][year]["granule_count"]
            lines.append(f"* NEW YEAR OF DATA: {year} added ({n} granule(s)).")

    if changes["removed_years"]:
        for year in changes["removed_years"]:
            lines.append(f"* Year {year} is no longer present in the active collection (was previously available).")

    for year, delta in sorted(changes["changed_years"].items()):
        parts = []
        if delta["added"]:
            parts.append(f"{len(delta['added'])} granule(s) added")
        if delta["removed"]:
            parts.append(f"{len(delta['removed'])} granule(s) removed")
        if delta["modified"]:
            parts.append(f"{len(delta['modified'])} granule(s) reprocessed/replaced (checksum, size, or revision changed)")
        lines.append(f"* Existing year {year} changed: {'; '.join(parts)}.")
        if delta["modified"]:
            shown = delta["modified"][:5]
            more = f" (+{len(delta['modified']) - 5} more)" if len(delta["modified"]) > 5 else ""
            lines.append(f"    e.g. {', '.join(shown)}{more}")

    return "\n".join(lines)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def check_dataset(short_name: str, state_file: Path | None = None, update: bool = False) -> tuple[str, dict]:
    """Fetch current CMR state for short_name and diff against state_file if it exists.

    Read-only by default: the state file is only created/overwritten when
    update=True is passed. Returns (human-readable report string, freshly fetched state dict).
    """
    state_file = state_file or default_state_file(short_name)

    new_state = build_state(short_name)

    old_state = None
    changes = None
    if state_file.exists():
        old_state = json.loads(state_file.read_text())
        changes = diff_states(old_state, new_state)

    report = format_report(old_state, new_state, changes, update)

    if update:
        state_file.write_text(json.dumps(new_state, indent=2, sort_keys=True))
        report += f"\n\n(state saved to {state_file})"
    else:
        report += "\n\n(no state file written — pass --update to save this state)"

    return report, new_state


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset", help="CMR short_name of the dataset to check, e.g. NSIDC-0731, NSIDC-0630, IDBMG4")
    parser.add_argument("--state-file", type=Path, default=None,
                         help="Where to cache state between runs (default: <dataset>_state.json next to this script)")
    parser.add_argument("--update", action="store_true",
                         help="Save the freshly fetched state to the state file (creating or overwriting it). "
                              "Default is read-only: report a diff (if a state file exists) but never write one.")
    parser.add_argument("--json", action="store_true", help="Also print the raw new-state JSON to stdout")
    args = parser.parse_args()

    try:
        state_file = args.state_file or default_state_file(args.dataset)
        report, new_state = check_dataset(args.dataset, state_file=state_file, update=args.update)
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
