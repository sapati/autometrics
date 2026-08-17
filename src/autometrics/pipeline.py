from astropy.wcs import WCS
from astropy.io import fits
from .calibration import build_masters, calibrate_lights
from .registration import align_lights, stack_images
from .astrometry import solve
from .catalog import fetch_sequence
from .detection import detect, estimate_fwhm
from .photometry import match_catalog, choose_comp_check, run_photometry, write_table
from .diagnostics import plot_field, plot_lightcurve
from .reporting import write_report


def run(cfg, target):
    print("Building calibration masters")
    darks, flats = build_masters(cfg)
    print("Calibrating lights")
    calibrated = calibrate_lights(cfg, darks, flats)
    print(f"Registering {len(calibrated)} frames")
    reference, aligned = align_lights(cfg, calibrated)
    stack = stack_images(cfg, aligned)
    print("Solving WCS")
    solved = solve(cfg, stack)
    wcs = WCS(fits.getheader(solved))
    detections = detect(solved, cfg)
    cfg.fwhm_pixels = estimate_fwhm(solved, detections, cfg)
    print(f"Detected {len(detections)} stars; estimated FWHM {cfg.fwhm_pixels:.2f} pixels")
    stars, chart_id = fetch_sequence(cfg, target)
    matches = match_catalog(wcs, stars, detections, cfg.match_radius_arcsec)
    if not any(m.star.role == "target" for m in matches):
        raise RuntimeError("Target was not matched to a detected star; check WCS and target visibility")
    comp, check = choose_comp_check(matches, aligned, cfg)
    rows = run_photometry(cfg, aligned, matches, comp, check)
    if not rows:
        raise RuntimeError("No usable photometry rows remained after SNR/flux checks")
    write_table(cfg.derived / "photometry" / "measurements.csv", rows)
    plot_field(cfg, solved, matches, comp, check)
    plot_lightcurve(cfg, rows)
    report = write_report(cfg, target, chart_id, rows, comp, check)
    print(f"Wrote {len(rows)} observations to {report}")
    return report
