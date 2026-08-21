# autometrics

`autometrics` is an aperture-photometry pipeline for CCD variable-star observations. Includes calibration, registration, local astap WCS solving, AAVSO VSP sequence matching, aperture selection, differential photometry, diagnostics, and AAVSO Extended Format export.

## Setup

```bash
conda create -n autometrics python=3.11 astropy photutils ccdproc astroalign matplotlib scipy requests typer
conda activate autometrics
python -m pip install -e .
cp config.sample.toml my-project/config.toml
```

Install `astap` and its index files separately. Configure `astrometry.solve_field` and, if needed, `index_path`.

## Layout

```text
my-project/source/{darks,flats,lights}
my-project/derived/{masterdark,masterflat,calibrated,astrometry,catalogs,photometry,figures,reports}
```

FITS headers should include `IMAGETYP`, `EXPTIME`, `FILTER`, and `DATE-OBS`. Observatory coordinates and observer code belong in `config.toml`.

## Run

```bash
autometrics my-project "SS Cyg" --config my-project/config.toml
```

The run writes an annotated PNG, machine-readable photometry tables, a JSON summary, and an AAVSO Extended Format CSV under `derived/`.
