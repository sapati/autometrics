import subprocess
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS


def solve(cfg, image):
    output = cfg.derived / "astrometry" / "solved.fits"
    if output.exists():
        return output
    command = [cfg.solve_field, "--overwrite", "--no-plots", "--new-fits", str(output), "--dir", str(cfg.derived / "astrometry")]
    if cfg.index_path:
        command += ["--index-dir", cfg.index_path]
    if cfg.pixel_scale:
        command += ["--scale-units", "arcsecperpix", "--scale-low", str(cfg.pixel_scale * .7), "--scale-high", str(cfg.pixel_scale * 1.3)]
    command.append(str(image))
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError("astrometry.net solve-field is not installed or not on PATH") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Plate solve failed: {exc.stderr[-1000:]}") from exc
    if not output.exists():
        candidates = list((cfg.derived / "astrometry").glob("*.new")) + list((cfg.derived / "astrometry").glob("*.fits"))
        if candidates:
            output = candidates[0]
        else:
            raise RuntimeError("solve-field completed without producing a solved FITS file")
    WCS(fits.getheader(output))
    return output
