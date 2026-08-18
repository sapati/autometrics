import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.optimize import least_squares


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
    """Estimate seeing from elliptical Moffat fits to detected stars."""
    data = fits.getdata(image).astype(float)
    values = []
    radius = max(4, int(3 * cfg.fwhm_pixels))
    for x, y in detections[:100]:
        xi, yi = int(round(x)), int(round(y))
        cut = data[max(0, yi-radius):yi+radius+1, max(0, xi-radius):xi+radius+1]
        if cut.size == 0 or not np.all(np.isfinite(cut)) or np.nanmax(cut) >= cfg.saturation:
            continue
        yy, xx = np.indices(cut.shape)
        background = np.median(cut)
        weights = np.clip(cut - background, 0, None)
        total = weights.sum()
        if total <= 0:
            continue
        x0 = np.sum(weights * xx) / total
        y0 = np.sum(weights * yy) / total
        amplitude = weights.max()
        scale = cfg.fwhm_pixels / (2 * np.sqrt(2 ** (1 / 3.5) - 1))
        span = max(np.ptp(cut), 1.0)

        def residuals(params):
            base, peak, xcen, ycen, alpha_x, alpha_y, beta, theta = params
            dx, dy = xx - xcen, yy - ycen
            xrot = np.cos(theta) * dx + np.sin(theta) * dy
            yrot = -np.sin(theta) * dx + np.cos(theta) * dy
            model = base + peak * (1 + (xrot / alpha_x) ** 2 + (yrot / alpha_y) ** 2) ** -beta
            return (model - cut).ravel()

        lower = [background - span, 0, 0, 0, 0.25, 0.25, 1.1, -np.pi / 2]
        upper = [np.nanmax(cut), np.inf, cut.shape[1] - 1, cut.shape[0] - 1,
                 max(cut.shape), max(cut.shape), 10, np.pi / 2]
        try:
            fit = least_squares(
                residuals,
                [background, amplitude, x0, y0, scale, scale, 3.5, 0],
                bounds=(lower, upper),
                loss="soft_l1",
                f_scale=max(np.std(cut - background), 1.0),
                max_nfev=400,
            )
        except (ValueError, RuntimeError):
            continue
        if not fit.success:
            continue
        alpha_x, alpha_y, beta = fit.x[4:7]
        factor = 2 * np.sqrt(2 ** (1 / beta) - 1)
        fwhm_x, fwhm_y = factor * alpha_x, factor * alpha_y
        fwhm = np.sqrt(fwhm_x * fwhm_y)
        if 1.0 <= fwhm_x <= 20.0 and 1.0 <= fwhm_y <= 20.0:
            values.append(fwhm)
    return float(np.median(values)) if values else cfg.fwhm_pixels
