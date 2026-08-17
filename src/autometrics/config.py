from dataclasses import dataclass
from pathlib import Path
import tomllib


@dataclass
class Config:
    root: Path
    observer_code: str = ""
    observer_name: str = ""
    software: str = "autometrics"
    latitude: float = 0.0
    longitude: float = 0.0
    elevation_m: float = 0.0
    gain: float = 1.0
    read_noise: float = 5.0
    saturation: float = 60000.0
    pixel_scale: float | None = None
    filter_name: str = "V"
    solve_field: str = "solve-field"
    index_path: str = ""
    radius_arcmin: float = 30.0
    detection_sigma: float = 5.0
    fwhm_pixels: float = 4.0
    aperture_fwhm: float = 1.8
    annulus_inner_fwhm: float = 3.0
    annulus_outer_fwhm: float = 5.0
    match_radius_arcsec: float = 3.0
    min_snr: float = 10.0
    dark_tolerance: float = 5.0

    @property
    def source(self): return self.root / "source"
    @property
    def darks(self): return self.source / "darks"
    @property
    def flats(self): return self.source / "flats"
    @property
    def lights(self): return self.source / "lights"
    @property
    def derived(self): return self.root / "derived"

    def ensure_dirs(self):
        for name in ("masterdark", "masterflat", "calibrated", "astrometry", "catalogs", "photometry", "figures", "reports"):
            (self.derived / name).mkdir(parents=True, exist_ok=True)
        (self.derived / "calibrated" / "aligned").mkdir(exist_ok=True)


def load_config(root: str | Path, filename: str | Path | None = None) -> Config:
    root = Path(root).expanduser().resolve()
    path = Path(filename) if filename else root / "config.toml"
    values = {}
    if path.is_file():
        with path.open("rb") as handle:
            raw = tomllib.load(handle)
        values = {
            "observer_code": raw.get("observer", {}).get("code", ""),
            "observer_name": raw.get("observer", {}).get("name", ""),
            "software": raw.get("observer", {}).get("software", "autometrics"),
            "latitude": raw.get("observer", {}).get("latitude", 0.0),
            "longitude": raw.get("observer", {}).get("longitude", 0.0),
            "elevation_m": raw.get("observer", {}).get("elevation_m", 0.0),
            "gain": raw.get("instrument", {}).get("gain_e_per_adu", 1.0),
            "read_noise": raw.get("instrument", {}).get("read_noise_e", 5.0),
            "saturation": raw.get("instrument", {}).get("saturation_adu", 60000.0),
            "pixel_scale": raw.get("instrument", {}).get("pixel_scale_arcsec"),
            "filter_name": raw.get("instrument", {}).get("filter", "V"),
            "solve_field": raw.get("astrometry", {}).get("solve_field", "solve-field"),
            "index_path": raw.get("astrometry", {}).get("index_path", ""),
            "radius_arcmin": raw.get("astrometry", {}).get("radius_arcmin", 30.0),
        }
        values.update({k: raw.get("photometry", {}).get(k, v) for k, v in {
            "detection_sigma": 5.0, "fwhm_pixels": 4.0, "aperture_fwhm": 1.8,
            "annulus_inner_fwhm": 3.0, "annulus_outer_fwhm": 5.0,
            "match_radius_arcsec": 3.0, "min_snr": 10.0}.items()})
        values["dark_tolerance"] = raw.get("processing", {}).get("dark_tolerance_seconds", 5.0)
    cfg = Config(root=root, **values)
    cfg.ensure_dirs()
    return cfg
