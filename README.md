# Green IBIS: Greenland Ice Bed from Ice Surface
Inferring seasonal timescale basal processes from ice surface observations and ice-flow models for Greenland outlet glaciers.

## Note: This repository is actively being built and will undergo many substantial changes.

## Repository Structure
This repository combines MATLAB and Python/Jupyter Notebook scripts to handle observations and models for all parts of the Green IBIS project, including the model-vs-observation comparison tools formerly in the ModelObsCompare (moc) repository. Code is organized by function, not by language: MATLAB and Python files that do the same job live side by side in the same `src/` subfolder.

The model-vs-observation tools (ROIs, flowlines, spatial/time-series comparison, variability maps) are documented in [docs/model_obs_comparison.md](docs/model_obs_comparison.md).

```
green-ibis/
├── setup_paths.m   # MATLAB: put src/ (and optionally ISSM) on the path
├── src/
│   ├── config/       # data-path resolution (MATLAB + config.py)
│   ├── physics/      # basal shear stress, driving stress, effective/water/ice pressure, ...
│   ├── conversions/  # mass <-> rate, volume -> mass, datenum <-> decimal year
│   ├── data_io/      # model output + observation readers/exporters
│   ├── roi/          # regions of interest: select, save, load, mask
│   ├── flowline/     # flowlines from shapefiles, clipping, sampling obs along them
│   ├── interp/       # model mesh -> grid / points / flowline
│   ├── compare/      # model vs obs: spatial differences, time series, variability
│   ├── plotting/     # maps, flowlines, time series, difference and variability plots
│   └── data_tools/   # dataset version-update checkers
├── workflows/      # one-off analysis drivers with hard-coded domains/paths
├── examples/       # portable demos of src/ functions (MATLAB + Python)
├── notebooks/      # Jupyter notebooks
├── tests/          # mirrors src/ subfolders
├── data/           # local static files, e.g. saved ROIs (not tracked by git)
├── docs/           # longer documentation
└── pyproject.toml  # Python dependencies (uv-managed)
```

A script goes in `workflows/` if you would need to edit values inside it (a domain name, catchment IDs, file paths) to reuse it. It goes in `src/` if it works for any domain or period through its inputs. Small, portable demonstrations of `src/` functions go in `examples/`.

## Getting started

### MATLAB
```matlab
cfg = setup_paths();          % adds src/, examples/, workflows/, tests/ (+ ISSM if configured)
configure_paths('-set')       % optional: prompt once for every data location
```

### Python
There is no installable package. Put `src/` on the import path and import modules by folder:
```bash
uv sync --extra shapefile --extra interactive   # create the environment from pyproject.toml / uv.lock
PYTHONPATH=src uv run python examples/example_sverdrups_workflow.py
```
```python
from data_io.obs_netcdf import load_obs_netcdf
from roi.roi import roi_from_obs
from compare.variance import velocity_variance
```
Notebooks and scripts in `notebooks/`, `examples/` and `tests/` add `src/` to `sys.path` themselves.

### Data paths
Nothing is hard-coded to one machine. MATLAB (`data_path`, `configure_paths`) and Python (`config.config`) share the same data-path settings. Each path is looked up in this order:

1. a session override (`configure_paths(name, value)` / `configure(name=value)`)
2. an environment variable: `ISSM_DIR`, `GREENIBIS_MODEL_DIR`, `GREENIBIS_OBS_FLOWLINE_DIR`, `GREENIBIS_OBS_NETCDF_DIR`, `GREENIBIS_FLOWLINE_SHP_DIR`, `GREENIBIS_TERMINI_SHP`, `GREENIBIS_ICEFRONTS_NC`, `GREENIBIS_ROI_DIR`
3. the shared config file `~/.config/green-ibis/paths.json` (override its location with `GREENIBIS_CONFIG`)
4. an interactive prompt, whose answer is saved to the config file

Show or set the paths from Python with `PYTHONPATH=src python -m config.config [--set]`.

### Tests
```matlab
smoke_test                        % needs configured data paths
test_variance_synthetic           % synthetic, no data needed
```
```bash
python tests/compare/test_variance_synthetic.py
python tests/flowline/test_flowline_shp_synthetic.py
python tests/data_io/test_model_synthetic.py
```
