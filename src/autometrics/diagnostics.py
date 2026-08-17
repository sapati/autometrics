import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, SqrtStretch
from photutils.aperture import CircularAperture, CircularAnnulus


def plot_field(cfg, image, matches, comp, check):
    data = fits.getdata(image)
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.imshow(data, origin="lower", cmap="gray", norm=ImageNormalize(data, stretch=SqrtStretch()))
    radius = cfg.aperture_fwhm * cfg.fwhm_pixels
    for match in matches:
        color = "red" if match.star.role == "target" else "lime" if match is comp else "cyan" if match is check else "yellow"
        CircularAperture((match.x, match.y), radius).plot(ax=ax, color=color, lw=1)
        CircularAnnulus((match.x, match.y), cfg.annulus_inner_fwhm * cfg.fwhm_pixels, cfg.annulus_outer_fwhm * cfg.fwhm_pixels).plot(ax=ax, color=color, lw=.7)
        ax.text(match.x + 4, match.y + 4, match.star.label, color=color, fontsize=9)
    ax.set_title("autometrics: target, comparison, check, and sequence stars")
    ax.set_xlabel("pixel x"); ax.set_ylabel("pixel y")
    output = cfg.derived / "figures" / "field.png"
    fig.savefig(output, dpi=720, bbox_inches="tight"); plt.close(fig)
    return output


def plot_lightcurve(cfg, rows):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar([r["jd"] for r in rows], [r["magnitude"] for r in rows], yerr=[r["error"] for r in rows], fmt=".")
    ax.invert_yaxis(); ax.set_xlabel("JD"); ax.set_ylabel("Magnitude")
    output = cfg.derived / "figures" / "lightcurve.png"
    fig.savefig(output, dpi=160, bbox_inches="tight"); plt.close(fig)
    return output
