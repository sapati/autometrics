from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats


def fits_files(path):
    return sorted(p for p in Path(path).glob("*.fit*") if p.is_file())


def _header(path):
    return fits.getheader(path)


def _type(header):
    return str(header.get("IMAGETYP", header.get("IMAGE-TYPE", ""))).lower()


def _exposure(header):
    return float(header.get("EXPTIME", header.get("EXPOSURE", 0.0)))


def _filter(header):
    return str(header.get("FILTER", "Unknown")).strip()


def _combine(paths):
    arrays = [fits.getdata(p).astype(float) for p in paths]
    return np.nanmedian(np.stack(arrays), axis=0)


def build_masters(cfg):
    darks = [p for p in fits_files(cfg.darks) if "dark" in _type(_header(p))]
    flats = [p for p in fits_files(cfg.flats) if "flat" in _type(_header(p))]
    dark_groups = {}
    for path in darks:
        dark_groups.setdefault(round(_exposure(_header(path)), 3), []).append(path)
    dark_map = {}
    for exposure, paths in dark_groups.items():
        output = cfg.derived / "masterdark" / f"combined_dark_{exposure:.3f}.fits"
        if not output.exists():
            data = _combine(paths)
            fits.writeto(output, data, _header(paths[0]), overwrite=True)
        dark_map[exposure] = output
    flat_groups = {}
    for path in flats:
        flat_groups.setdefault(_filter(_header(path)), []).append(path)
    flat_map = {}
    for filt, paths in flat_groups.items():
        output = cfg.derived / "masterflat" / f"combined_flat_filter_{filt}.fits"
        if not output.exists():
            corrected = []
            for path in paths:
                data = fits.getdata(path).astype(float)
                hdr = _header(path)
                if dark_map:
                    dark = min(dark_map, key=lambda x: abs(x - _exposure(hdr)))
                    data -= fits.getdata(dark)
                _, med, _ = sigma_clipped_stats(data)
                corrected.append(data / med if med else data)
            fits.writeto(output, np.nanmedian(corrected, axis=0), _header(paths[0]), overwrite=True)
        flat_map[filt] = output
    return dark_map, flat_map


def calibrate_lights(cfg, dark_map, flat_map):
    outputs = []
    for path in fits_files(cfg.lights):
        hdr = _header(path)
        if "light" not in _type(hdr):
            continue
        output = cfg.derived / "calibrated" / path.name
        if not output.exists():
            data = fits.getdata(path).astype(float)
            if dark_map:
                dark = min(dark_map, key=lambda x: abs(x - _exposure(hdr)))
                if abs(dark - _exposure(hdr)) > cfg.dark_tolerance:
                    raise RuntimeError(f"No dark within tolerance for {path.name}")
                data -= fits.getdata(dark)
            flat = flat_map.get(_filter(hdr))
            if flat:
                data /= np.where(fits.getdata(flat) == 0, 1, fits.getdata(flat))
            hdr["AUTOCAL"] = (True, "Calibrated by autometrics")
            fits.writeto(output, data, hdr, overwrite=True)
        outputs.append(output)
    return outputs
