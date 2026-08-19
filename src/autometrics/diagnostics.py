import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, PercentileInterval, SqrtStretch
from photutils.aperture import CircularAperture, CircularAnnulus


def plot_field(cfg, image, matches, comp, check):
    data = fits.getdata(image)
    radius = cfg.aperture_fwhm * cfg.fwhm_pixels
    relevant = [match for match in matches if match.star.role == "target" or match is comp or match is check]
    xs = [match.x for match in relevant]
    ys = [match.y for match in relevant]
    margin = max(4 * radius, 0.1 * max(max(xs) - min(xs), max(ys) - min(ys)))
    x_min = max(0, min(xs) - margin)
    x_max = min(data.shape[1], max(xs) + margin)
    y_min = max(0, min(ys) - margin)
    y_max = min(data.shape[0], max(ys) + margin)

    width, height = x_max - x_min, y_max - y_min
    fig, ax = plt.subplots(figsize=(12, 12 * height / width))
    ax.set_position([0, 0, 1, 1])
    ax.imshow(
        data,
        origin="lower",
        cmap="gray_r",
        norm=ImageNormalize(data, interval=PercentileInterval(99.5), stretch=SqrtStretch()),
    )
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_axis_off()
    for match in relevant:
        color = "red" if match.star.role == "target" else "lime" if match is comp else "cyan" if match is check else "yellow"
        CircularAperture((match.x, match.y), radius).plot(ax=ax, color=color, lw=1)
        CircularAnnulus((match.x, match.y), cfg.annulus_inner_fwhm * cfg.fwhm_pixels, cfg.annulus_outer_fwhm * cfg.fwhm_pixels).plot(ax=ax, color=color, lw=.7)
        ax.text(match.x + 4, match.y + 4, match.star.label, color=color, fontsize=9)
    output = cfg.derived / "figures" / "field.png"
    fig.savefig(output, dpi=150, pad_inches=0)
    plt.close(fig)
    return output


def plot_lightcurve(cfg, rows):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar([r["jd"] for r in rows], [r["magnitude"] for r in rows], yerr=[r["error"] for r in rows], fmt=".")
    ax.invert_yaxis(); ax.set_xlabel("JD"); ax.set_ylabel("Magnitude")
    output = cfg.derived / "figures" / "lightcurve.png"
    fig.savefig(output, dpi=160, bbox_inches="tight"); plt.close(fig)
    return output
