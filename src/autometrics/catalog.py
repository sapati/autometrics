import json
from pathlib import Path
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from .models import SequenceStar


def _deg(value, ra=False):
    return SkyCoord(value, unit=(u.hourangle if ra else u.deg)).degree


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
    target_star = SequenceStar(payload.get("auid", target), target, _deg(payload["ra"], True), _deg(payload["dec"]), None, None, role="target")
    stars = [target_star]
    for item in payload.get("photometry", []):
        bands = {b["band"].upper(): b for b in item.get("bands", [])}
        band = bands.get(cfg.filter_name.upper()) or bands.get("V")
        if not band or band.get("mag") is None:
            continue
        stars.append(SequenceStar(item["auid"], item.get("label", item["auid"]), _deg(item["ra"], True), _deg(item["dec"]), float(band["mag"]), float(band.get("error") or 0), item.get("comments") or ""))
    return stars, payload.get("chartid", "")
