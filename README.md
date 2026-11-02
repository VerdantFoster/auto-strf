# STRF Automation

STRF Automation processes SatNOGS radio observations and generates candidate Two-Line Element (TLE) updates using the orbit-fitting engine from [STRF](https://github.com/cbassa/strf). It downloads waterfall images and metadata, extracts signal tracks, reconstructs Doppler measurements, and produces fitted TLEs, residuals, and plots.

## Installation

On Debian or Ubuntu, install the build dependencies:

```sh
sudo apt update
sudo apt install make gcc gfortran pgplot5 libpng-dev libx11-dev libgsl-dev libfftw3-dev libsox-dev python3-venv
```

If PGPLOT headers are missing, install your distribution's PGPLOT development package.

From the repository directory, create a Python environment and build the tools:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r satnogs_waterfall/requirements.txt
make
make rffit-batch
```

Run commands from this directory. Activate the environment again with `. .venv/bin/activate` when opening a new terminal. Downloads require internet access; interactive marking requires a graphical desktop.

## Usage

### Process an observation

Replace the example ID with a completed SatNOGS observation that has a waterfall image and tracking TLE:

```sh
OBS_ID=12345678
python -m satnogs_waterfall fetch "$OBS_ID"
python -m satnogs_waterfall inspect "$OBS_ID"
python -m satnogs_waterfall suggest "$OBS_ID"
python -m satnogs_waterfall export "$OBS_ID" --source proposal
python -m satnogs_waterfall.fit_outputs "$OBS_ID" --mode mamm
```

The commands download the observation, inspect its calibration, propose signal points, export Doppler measurements, and fit a candidate TLE. Station IDs must be between 1 and 9999. Automatically extracted points remain unreviewed until inspected by the user.

### Fetch the latest observation

Specify a NORAD catalog number and one or more SatNOGS station IDs:

```sh
NORAD_ID=25544
STATION_ID=1
python -m satnogs_waterfall fetch-latest --norad-id "$NORAD_ID" --stations "$STATION_ID"
```

Replace these example values with your target satellite and station. Use the observation ID printed by the command as `OBS_ID`, then continue from `inspect` in the workflow above.

### Review signal points

To review or adjust the proposed points before repeating the fit:

```sh
python -m satnogs_waterfall mark "$OBS_ID"
python -m satnogs_waterfall export "$OBS_ID" --source selection
python -m satnogs_waterfall.fit_outputs "$OBS_ID" --mode mamm --replace
```

In the marking window, click to add a point, press `u` to undo, `c` to clear, and `f` to save. Use `--replace` when regenerating existing fit outputs.

The default frequency reconstruction assumes a Doppler-corrected SatNOGS waterfall and uses its tracking TLE. For a waterfall known to be uncorrected, specify the correction state when exporting:

```sh
python -m satnogs_waterfall export "$OBS_ID" --source selection \
  --correction uncorrected --evidence "Receiver recorded without Doppler correction"
```

### Choose a fitting mode

| Mode | Adjusted orbital elements |
| --- | --- |
| `fixed` | None |
| `ma` | Mean anomaly |
| `mamm` | Mean anomaly and mean motion |
| `raan` | Right ascension of the ascending node (RAAN) |
| `maraan` | Mean anomaly and RAAN |
| `mammraan` | Mean anomaly, mean motion, and RAAN |

The transmitter frequency is estimated in every mode, including `fixed`.

### Read the outputs

Files are saved in `data/waterfall/<observation_id>/`:

- `measurements.dat`: exported Doppler measurements for STRF.
- `seed.tle`: tracking TLE used to initialize the fit.
- `rffit_candidate.tle`: fitted candidate TLE.
- `rffit_residuals.csv`: observed frequencies, predictions, and residuals.
- `rffit_result.png`: frequency fit and residual panels on a white background.
- `rffit_full.png`: complete STRF view on a black background, including the Doppler curve, radial-velocity axis, sky plot, candidate TLE, measurement count, RMS, TCA (when available), TLE epoch, and fitted-parameter indicators.
- `rffit_provenance.json`: input hashes, fitting settings, and metrics.
- `rffit.log`: fitting log.

Both PNGs are saved automatically by `satnogs_waterfall.fit_outputs`. No graphical desktop is needed for file export. The full view uses the saved candidate TLE and the transmitter frequency from the batch fit; it does not run another optimization. Parameter indicators show the selected fitting mode and are not interactive in the saved image. The full-view RMS is checked against the residual CSV before the outputs are saved.

To update an existing installation and regenerate the plots, run from the repository directory with your Python environment activated:

```sh
git pull --ff-only
make rffit rffit-batch
python -m satnogs_waterfall.fit_outputs "$OBS_ID" --mode mamm --replace
```

Set `OBS_ID` to an observation that has already been exported, and choose the same fitting mode you intend to use. `--replace` repeats the fit and replaces its existing outputs, including the candidate TLE and both plots. Successful waterfall calibration and measurement export are still required.

Candidate TLEs do not automatically replace an operational catalog. Review the signal selection and validate predictions against independent observations before using a candidate for tracking.

For command options or a command to open the measurements in interactive STRF:

```sh
python -m satnogs_waterfall --help
python -m satnogs_waterfall.fit_outputs --help
python -m satnogs_waterfall fit-command "$OBS_ID"
```
