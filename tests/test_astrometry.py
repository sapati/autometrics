import numpy as np
from astropy.io import fits

from autometrics.astrometry import _has_data, _merge_image_and_wcs


def test_wcs_header_is_merged_with_original_pixels(tmp_path):
    image = tmp_path / "image.fits"
    solved = tmp_path / "solved.fits"
    data = np.arange(16, dtype=np.float32).reshape(4, 4)
    fits.writeto(image, data, fits.Header({"FILTER": "V"}))
    solved_header = fits.Header({"CRVAL1": 10.0, "CRVAL2": 20.0, "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"})

    _merge_image_and_wcs(image, solved_header, solved)

    assert _has_data(solved)
    with fits.open(solved) as hdul:
        assert np.array_equal(hdul[0].data, data)
        assert hdul[0].header["CRVAL1"] == 10.0
        assert hdul[0].header["FILTER"] == "V"
