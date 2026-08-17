import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder


def detect(image, cfg):
    data = fits.getdata(image).astype(float)
    _, median, std = sigma_clipped_stats(data)
    finder = DAOStarFinder(fwhm=cfg.fwhm_pixels, threshold=cfg.detection_sigma * std)
    table = finder(data - median)
    if table is None:
        return np.empty((0, 2))
    pad = int(3 * cfg.fwhm_pixels)
    good = ((table["xcentroid"] > pad) & (table["xcentroid"] < data.shape[1] - pad) &
            (table["ycentroid"] > pad) & (table["ycentroid"] < data.shape[0] - pad) &
            (table["peak"] < cfg.saturation))
    return np.column_stack((np.asarray(table["xcentroid"])[good], np.asarray(table["ycentroid"])[good]))


def estimate_fwhm(image, detections, cfg):
    """Estimate seeing from weighted second moments of isolated detections."""
    data = fits.getdata(image).astype(float)
    values = []
    radius = max(3, int(2.5 * cfg.fwhm_pixels))
    for x, y in detections[:100]:
        xi, yi = int(round(x)), int(round(y))
        cut = data[max(0, yi-radius):yi+radius+1, max(0, xi-radius):xi+radius+1]
        if cut.size == 0:
            continue
        yy, xx = np.indices(cut.shape)
        weights = np.clip(cut - np.nanmedian(cut), 0, None)
        total = weights.sum()
        if total <= 0:
            continue
        sx = np.sqrt(np.sum(weights * (xx - np.sum(weights * xx) / total) ** 2) / total)
        sy = np.sqrt(np.sum(weights * (yy - np.sum(weights * yy) / total) ** 2) / total)
        fwhm = 2.355 * np.sqrt((sx * sx + sy * sy) / 2)
        if 1.0 <= fwhm <= 20.0:
            values.append(fwhm)
    return float(np.median(values)) if values else cfg.fwhm_pixels
