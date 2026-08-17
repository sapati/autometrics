from dataclasses import dataclass


@dataclass
class SequenceStar:
    auid: str
    label: str
    ra_deg: float
    dec_deg: float
    magnitude: float | None
    mag_error: float | None
    comments: str = ""
    role: str = "sequence"


@dataclass
class Match:
    star: SequenceStar
    x: float
    y: float
    distance_arcsec: float
    fwhm: float
    isolated: bool = True
