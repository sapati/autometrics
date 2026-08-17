from pathlib import Path
import csv
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from autometrics.config import load_config
from autometrics.models import SequenceStar
from autometrics.photometry import match_catalog, write_table
from autometrics.reporting import write_report
from autometrics.calibration import calibrate_lights
from autometrics.catalog import _coords


def test_config_creates_reference_layout(tmp_path):
    cfg = load_config(tmp_path)
    assert cfg.darks == tmp_path / "source" / "darks"
    assert (tmp_path / "derived" / "reports").is_dir()


def test_wcs_catalog_match():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [50, 50]
    wcs.wcs.cdelt = [-1 / 3600, 1 / 3600]
    wcs.wcs.crval = [10, 20]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    star = SequenceStar("A", "A", 10, 20, 12.0, 0.02, role="target")
    matches = match_catalog(wcs, [star], np.array([[49, 49]]), 3)
    assert len(matches) == 1
    assert matches[0].star.auid == "A"


def test_report_has_required_header(tmp_path):
    cfg = load_config(tmp_path)
    comp = type("M", (), {"star": type("S", (), {"auid": "C", "magnitude": 12.3})()})()
    check = type("M", (), {"star": type("S", (), {"auid": "K", "magnitude": 13.1})()})()
    rows = [{"jd": 2450000.1, "magnitude": 12.1, "error": .02, "target_flux": 1000, "comp_flux": 2000, "check_magnitude": 13.0, "airmass": 1.2}]
    output = write_report(cfg, "Target", "X1", rows, comp, check)
    text = output.read_text()
    assert "#TYPE=Extended" in text
    assert "#OBSCODE=" in text
    assert "Target,2450000.10000" in text


def test_calibration_resolves_dark_exposure_key(tmp_path):
    cfg = load_config(tmp_path)
    dark_path = cfg.derived / "masterdark" / "dark.fits"
    light_path = cfg.lights / "light.fits"
    header = fits.Header({"IMAGETYP": "Light Frame", "EXPTIME": 10.0, "FILTER": "V"})
    fits.writeto(dark_path, np.ones((4, 4)), header, overwrite=True)
    fits.writeto(light_path, np.full((4, 4), 5.0), header, overwrite=True)
    output = calibrate_lights(cfg, {10.0: dark_path}, {})
    assert len(output) == 1
    assert np.allclose(fits.getdata(output[0]), 4.0)


def test_aavso_sexagesimal_coordinates():
    ra, dec = _coords("22:01:42.86", "69:44:36.5")
    assert 330 < ra < 331
    assert 69 < dec < 70
