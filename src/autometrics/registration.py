from pathlib import Path
import shutil
import numpy as np
import astroalign as aa
from astropy.io import fits


def align_lights(cfg, paths):
    if not paths:
        raise RuntimeError("No calibrated light frames found")
    reference = paths[0]
    aligned = []
    ref_data = fits.getdata(reference).astype(float)
    ref_header = fits.getheader(reference)
    for path in paths:
        output = cfg.derived / "calibrated" / "aligned" / path.name
        if not output.exists():
            if path == reference:
                shutil.copy2(path, output)
            else:
                source_data = np.asarray(fits.getdata(path), dtype=np.float64)
                transform, _ = aa.find_transform(source_data, ref_data)
                data, footprint = aa.apply_transform(transform, source_data, ref_data)
                hdr = fits.getheader(path)
                for key in ("WCSAXES", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"):
                    if key in ref_header:
                        hdr[key] = ref_header[key]
                hdr["AUTOALGN"] = (True, "Registered by autometrics")
                fits.writeto(output, np.asarray(data, dtype=float), hdr, overwrite=True)
        aligned.append(output)
    return reference, aligned


def stack_images(cfg, paths):
    output = cfg.derived / "calibrated" / "stack.fits"
    if not output.exists():
        data = np.nanmedian(np.stack([fits.getdata(p).astype(float) for p in paths]), axis=0)
        fits.writeto(output, data, fits.getheader(paths[0]), overwrite=True)
    return output
