import subprocess
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS


def solve(cfg, image):
    output = cfg.derived / "astrometry" / "solved"
    output_fn = cfg.derived / "astrometry" / "solved.wcs"
    if output_fn.exists():
        return image
    command = [cfg.solve_field, "-f", str(image), "-o", str(output), "-update" ]
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
    WCS(fits.getheader(image))
    return image
