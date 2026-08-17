import json
from pathlib import Path
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from .models import SequenceStar


def _coords(ra, dec):
    """Parse AAVSO decimal or sexagesimal coordinates as a coordinate pair."""
    if not str(ra).strip() or not str(dec).strip():
        raise ValueError(f"AAVSO returned empty coordinates: ra={ra!r}, dec={dec!r}")
    coordinate = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
    return coordinate.ra.degree, coordinate.dec.degree


def fetch_sequence(cfg, target):
    cache = cfg.derived / "catalogs" / (target.replace("/", "_").replace(" ", "_") + ".json")
    if cache.exists():
        payload = json.loads(cache.read_text())
    else:
        url = "https://app.aavso.org/vsp/api/chart/"
        response = requests.get(url, params={"star": target, "fov": cfg.radius_arcmin, "maglimit": 18, "format": "json"}, timeout=30)
        response.raise_for_status()
        payload = response.json()
        cache.write_text(json.dumps(payload, indent=2))
    target_ra, target_dec = _coords(payload.get("ra", ""), payload.get("dec", ""))
    target_star = SequenceStar(payload.get("auid", target), target, target_ra, target_dec, None, None, role="target")
    stars = [target_star]
    for item in payload.get("photometry", []):
        bands = {b["band"].upper(): b for b in item.get("bands", [])}
        band = bands.get(cfg.filter_name.upper()) or bands.get("V")
        if not band or band.get("mag") is None:
            continue
        try:
            ra, dec = _coords(item.get("ra", ""), item.get("dec", ""))
        except ValueError:
            continue
        stars.append(SequenceStar(item["auid"], item.get("label", item["auid"]), ra, dec, float(band["mag"]), float(band.get("error") or 0), item.get("comments") or ""))
    return stars, payload.get("chartid", "")
