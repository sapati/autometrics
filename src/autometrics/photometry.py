import csv
import math
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from .models import Match


def match_catalog(wcs, stars, detections, radius_arcsec):
    if not len(detections):
        raise RuntimeError("No stars detected")
    world = wcs.pixel_to_world(detections[:, 0], detections[:, 1])
    # WCS can return a tuple for non-celestial or multi-coordinate WCS objects.
    if isinstance(world, tuple):
        if len(world) < 2:
            raise RuntimeError("Solved WCS did not provide celestial coordinates")
        world = SkyCoord(world[0], world[1], unit=(u.deg, u.deg), frame="icrs")
    elif not isinstance(world, SkyCoord):
        world = SkyCoord(world.ra, world.dec, frame="icrs")
    result = []
    for star in stars:
        target = SkyCoord(ra=star.ra_deg * u.deg, dec=star.dec_deg * u.deg, frame="icrs")
        distances = target.separation(world).arcsec
        index = int(np.argmin(distances))
        if distances[index] <= radius_arcsec:
            result.append(Match(star, float(detections[index, 0]), float(detections[index, 1]), float(distances[index]), 0.0))
    return result


def choose_comp_check(matches, frames, cfg):
    candidates = [m for m in matches if m.star.role == "sequence" and m.star.magnitude is not None and "do not use" not in m.star.comments.lower()]
    if len(candidates) < 2:
        raise RuntimeError("At least two usable AAVSO sequence stars are required for comparison and check")
    scores = []
    for match in candidates:
        values = []
        for frame in frames:
            values.append(measure(frame, match.x, match.y, cfg)[0])
        values = np.asarray(values)
        stable = np.nanstd(values) if np.isfinite(values).any() else 1e9
        snr_penalty = 0 if np.nanmedian(values) >= cfg.min_snr else 100
        scores.append((stable + snr_penalty + (match.star.mag_error or 0), match))
    scores.sort(key=lambda item: item[0])
    return scores[0][1], scores[1][1]


def measure(image, x, y, cfg):
    data = fits.getdata(image).astype(float)
    r = cfg.aperture_fwhm * cfg.fwhm_pixels
    rin = cfg.annulus_inner_fwhm * cfg.fwhm_pixels
    rout = cfg.annulus_outer_fwhm * cfg.fwhm_pixels
    position = [(x, y)]
    aperture = CircularAperture(position, r=r)
    annulus = CircularAnnulus(position, r_in=rin, r_out=rout)
    source_result = aperture_photometry(data, aperture)
    sky_result = aperture_photometry(data, annulus)
    source_sum = float(source_result["aperture_sum"][0])
    sky_sum = float(sky_result["aperture_sum"][0])
    sky_per_pixel = sky_sum / max(annulus.area, 1)
    flux = source_sum - sky_per_pixel * aperture.area
    source_variance = max(source_sum, 0) * cfg.gain
    sky_variance = max(sky_sum, 0) * (aperture.area / max(annulus.area, 1))
    error = math.sqrt(source_variance + sky_variance + aperture.area * (cfg.read_noise ** 2)) / cfg.gain
    snr = flux / error if error > 0 else 0
    return flux, error, snr, r, rin, rout, sky_per_pixel


def run_photometry(cfg, frames, matches, comp, check):
    target = next(m for m in matches if m.star.role == "target")
    rows = []
    for frame in frames:
        tflux, terr, tsnr, *_ = measure(frame, target.x, target.y, cfg)
        cflux, cerr, _, *_ = measure(frame, comp.x, comp.y, cfg)
        kflux, kerr, _, *_ = measure(frame, check.x, check.y, cfg)
        if min(tflux, cflux, kflux) <= 0 or min(tsnr, cfg.min_snr) < cfg.min_snr:
            continue
        differential = -2.5 * math.log10(tflux / cflux)
        magnitude = comp.star.magnitude + differential
        error = 1.0857 * math.sqrt((terr / tflux) ** 2 + (cerr / cflux) ** 2 + (comp.star.mag_error or 0) ** 2)
        check_mag = comp.star.magnitude + (-2.5 * math.log10(kflux / cflux))
        header = fits.getheader(frame)
        date = header.get("DATE-OBS")
        if not date:
            raise RuntimeError(f"DATE-OBS is missing from {frame.name}")
        from astropy.time import Time
        jd = Time(date, format="isot", scale="utc").jd
        rows.append({"file": frame.name, "jd": jd, "magnitude": magnitude, "error": error, "check_magnitude": check_mag, "comp_flux": cflux, "target_flux": tflux, "check_flux": kflux, "airmass": header.get("AIRMASS", "na")})
    return rows


def write_table(path, rows):
    fields = list(rows[0]) if rows else ["file", "jd", "magnitude", "error"]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
