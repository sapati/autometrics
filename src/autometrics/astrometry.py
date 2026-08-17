import subprocess
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS


def _merge_image_and_wcs(image, solved_header, output):
    """Write the original pixels with the astrometry.net WCS header."""
    with fits.open(image, memmap=False) as original:
        data = original[0].data
        if data is None:
            raise RuntimeError(f"Original image has no primary image data: {image}")
        header = original[0].header.copy()
    header.update(solved_header)
    header["AUTOWCS"] = (True, "WCS merged by autometrics")
    fits.PrimaryHDU(data=data, header=header).writeto(output, overwrite=True)


def _has_data(path):
    try:
        with fits.open(path, memmap=False) as hdul:
            return hdul[0].data is not None
    except (OSError, IndexError):
        return False


def _read_solved_header(path):
    """Read either astrometry.net's text WCS header or a FITS HDU."""
    if str(path).endswith(".wcs"):
        return fits.Header.fromtextfile(path)
    try:
        with fits.open(path, memmap=False) as solved:
            return solved[0].header.copy()
    except OSError:
        return fits.Header.fromtextfile(path)


def solve(cfg, image):
    output = cfg.derived / "astrometry" / "solved"
    output_fn = cfg.derived / "astrometry" / "solved.wcs"
    merged = cfg.derived / "astrometry" / "solved.fits"
    if merged.exists() and _has_data(merged):
        return merged
    if not output_fn.exists():
        command = [cfg.solve_field, "-f", str(image), "-o", str(output) ]
        try:
            subprocess.run(command, check=True, capture_output=True, text=True, timeout=cfg.solve_timeout_seconds)
        except FileNotFoundError as exc:
            raise RuntimeError("astrometry.net solve-field is not installed or not on PATH") from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"Plate solve failed: {exc.stderr[-1000:]}") from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(f"Plate solve timed out after {cfg.solve_timeout_seconds:.0f} seconds") from exc
    if not output_fn.exists():
        raise RuntimeError("solve-field completed without producing a solved FITS file")
    solved_header = _read_solved_header(output_fn)
    _merge_image_and_wcs(image, solved_header, merged)
    WCS(fits.getheader(merged))
    return merged
