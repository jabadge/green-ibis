#!/usr/bin/env python3
"""
check_greenibis_dataset_updates.py

Wrapper that checks external observational datasets referenced by
runme_greenibis.m for updates, by driving the two dataset-specific checkers
in this same directory:

    check_cmr_dataset_updates.py       (NASA CMR datasets)
    check_racmo_fgrn055_updates.py     (RACMO2.3p2 FGRN055)

Usage:
    python3 check_cryo23_dataset_updates.py                # report only, no state files written
    python3 check_cryo23_dataset_updates.py --update       # report, then save/update every state file
    python3 check_cryo23_dataset_updates.py --json         # also print each dataset's raw state as JSON

No third-party dependencies required (both wrapped scripts use urllib from
the standard library; this wrapper imports them directly).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SOURCE_SCRIPT = "runme_greenibis.m"


def _import_from_path(module_name: str, file_name: str):
    path = SCRIPT_DIR / file_name
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cmr_checker = _import_from_path("check_cmr_dataset_updates", "check_cmr_dataset_updates.py")
racmo_checker = _import_from_path("check_racmo_fgrn055_updates", "check_racmo_fgrn055_updates.py")


# --------------------------------------------------------------------------
# Dataset registry — see module docstring for how this was derived
# --------------------------------------------------------------------------

CMR_DATASETS = [
    {
        "short_name": "IDBMG4",
        "issm_function": "interpBedmachineGreenland.m",
        "note": "BedMachine Greenland v6.6",
    },
    {
        "short_name": "NSIDC-0645",
        "issm_function": "interpGimpdem.m",
        "note": "GIMP 90 m DEM (reads 'gimpdem_90m.tif')",
    },
    {
        "short_name": "NSIDC-0670",
        "issm_function": "interpJoughinCompositeGreenland.m; interpFromMEaSUREsGeotiffs.m (product 670)",
        "note": "MEaSUREs GrIMP Multi-year Greenland Ice Sheet Velocity Mosaic",
    },
    {
        "short_name": "NSIDC-0481",
        "issm_function": "interpFromMEaSUREsGeotiffs.m (product 481)",
        "note": "MEaSUREs GrIMP",
    },
    {
        "short_name": "NSIDC-0731",
        "issm_function": "interpFromMEaSUREsGeotiffs.m (product 731)",
        "note": "MEaSUREs GrIMP",
    },
    {
        "short_name": "NSIDC-0793",
        "issm_function": "interpMonthlyIceMaskGreene.m",
        "note": "greenland-icemask Greene et al. (2024)",
    },
]

IGNORED_ITEMS = [
    {
        "item": "ISMIP6Greenland_Thermal.mat",
        "used_by": "runme_greenibis.m directly, via loadmodel() (thermal spin-up state)",
        "reason": "ISMIP6 initMIP model output, distributed through the ISMIP6 project's own data portal.",
    },
    {
        "item": "ISMIP6 atmosphere forcing (aSMB_observed, e.g. MIROC5-rcp26)",
        "used_by": "interpISMIP6GreenlandSMB.m",
        "reason": "Distributed via the ISMIP6 project data portal (Globus/OSF).",
    },
    {
        "item": "ISMIP6 retreat masks",
        "used_by": "interpISMIP6GreenlandRetreatMask.m",
        "reason": "Same ISMIP6 project data portal as the atmosphere forcing above.",
    },
    {
        "item": "SeaRISE Greenland_5km geothermal heat flux (Shapiro-Ritzwoller)",
        "used_by": "interpSeaRISE.m",
        "reason": "Hosted at ciei.colorado.edu (University of Colorado).",
    },
    {
        "item": "RACMO2 1 km SMB_MEAN1960-1989_150m.nc (static 1960-1989 mean climatology)",
        "used_by": "interpRACMO1km.m",
        "reason": "A single static climatology file, not the monthly RACMO2.3p2 FGRN055 time series "
                  "check_racmo_fgrn055_updates.py tracks, and not clearly identifiable as the IMAU "
                  "page's '1 km downscaled' product either.",
    },
]


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def run(update: bool, min_resolution_km: float, state_dir: Path | None, want_json: bool) -> int:
    print(f"Dataset update check for {SOURCE_SCRIPT}")
    print("=" * (24 + len(SOURCE_SCRIPT)))
    print()

    all_states = {}
    errors = []

    print(f"--- NASA CMR datasets ({len(CMR_DATASETS)}) ---\n")
    for entry in CMR_DATASETS:
        short_name = entry["short_name"]
        state_file = (state_dir / f"{short_name}_state.json") if state_dir else cmr_checker.default_state_file(short_name)
        print(f"[{short_name}] used by {entry['issm_function']}")
        print(f"  {entry['note']}")
        try:
            report, state = cmr_checker.check_dataset(short_name, state_file=state_file, update=update)
            all_states[short_name] = state
        except RuntimeError as e:
            print(f"  ERROR: {e}")
            errors.append((short_name, str(e)))
            print()
            continue
        for line in report.splitlines():
            print(f"  {line}")
        print()

    print(f"--- RACMO2.3p2 FGRN055 ---\n")
    print("  used by interpRACMO23p2smb.m")
    racmo_state_file = (state_dir / "racmo_fgrn055_state.json") if state_dir else racmo_checker.DEFAULT_STATE_FILE
    try:
        report, state = racmo_checker.check_racmo(
            racmo_checker.DEFAULT_ZENODO_RECORD, racmo_checker.UU_PAGE_URL, racmo_state_file,
            min_resolution_km, update,
        )
        all_states["RACMO2.3p2_FGRN055"] = state
        for line in report.splitlines():
            print(f"  {line}")
    except RuntimeError as e:
        print(f"  ERROR: {e}")
        errors.append(("RACMO2.3p2_FGRN055", str(e)))
    print()

    print(f"--- Ignored — referenced by runme_Cryo23.m but not in CMR or the RACMO Zenodo archive ({len(IGNORED_ITEMS)}) ---")
    print("(listed so you can double-check this determination)\n")
    for entry in IGNORED_ITEMS:
        print(f"* {entry['item']}")
        print(f"    used by: {entry['used_by']}")
        print(f"    reason:  {entry['reason']}")
    print()

    if errors:
        print(f"--- {len(errors)} dataset(s) failed to check ---")
        for name, msg in errors:
            print(f"* {name}: {msg}")
        print()

    if want_json:
        print("--- raw state (JSON) ---")
        print(json.dumps(all_states, indent=2, sort_keys=True))

    return 1 if errors else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--update", action="store_true",
                         help="Save the freshly fetched state for every dataset (creating or overwriting each "
                              "state file). Default is read-only: report diffs where a state file already "
                              "exists, but never write one.")
    parser.add_argument("--min-resolution-km", type=float, default=racmo_checker.DEFAULT_MIN_RESOLUTION_KM,
                         help=f"Passed through to the RACMO checker (default: {racmo_checker.DEFAULT_MIN_RESOLUTION_KM})")
    parser.add_argument("--state-dir", type=Path, default=None,
                         help="Directory for all state files (default: next to each checker script, i.e. this "
                              "script's directory)")
    parser.add_argument("--json", action="store_true", help="Also print every dataset's raw fetched state as JSON")
    args = parser.parse_args()

    if args.state_dir:
        args.state_dir.mkdir(parents=True, exist_ok=True)

    return run(args.update, args.min_resolution_km, args.state_dir, args.json)


if __name__ == "__main__":
    raise SystemExit(main())
