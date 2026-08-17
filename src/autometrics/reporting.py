import csv
import json
import math
from pathlib import Path


def write_report(cfg, target, chart_id, rows, comp, check):
    output = cfg.derived / "reports" / f"{target.replace(' ', '_')}_aavso.csv"
    columns = ["NAME", "DATE", "MAG", "MERR", "FILT", "TRANS", "MTYPE", "CNAME", "CMAG", "KNAME", "KMAG", "AMASS", "GROUP", "CHART", "NOTES"]
    with output.open("w", newline="") as handle:
        handle.write("#TYPE=Extended\n#OBSCODE={0}\n#SOFTWARE={1}\n#DELIM=,\n#DATE=JD\n#OBSTYPE=CCD\n".format(cfg.observer_code, cfg.software))
        handle.write("#" + ",".join(columns) + "\n")
        writer = csv.writer(handle)
        for row in rows:
            comp_instrumental = -2.5 * math.log10(row["comp_flux"])
            check_instrumental = -2.5 * math.log10(row["check_flux"])
            notes = f"|VMAGINS={-2.5 * math.log10(row['target_flux']):.5f}|CMAGINS={comp_instrumental:.5f}|KMAGINS={check_instrumental:.5f}|CREFMAG={comp.star.magnitude}"
            writer.writerow([target, f"{row['jd']:.5f}", f"{row['magnitude']:.3f}", f"{row['error']:.3f}", cfg.filter_name, "NO", "STD", comp.star.auid, f"{comp_instrumental:.3f}", check.star.auid, f"{check_instrumental:.3f}", row["airmass"], "na", chart_id or "na", notes])
    summary = cfg.derived / "reports" / "run_summary.json"
    summary.write_text(json.dumps({"target": target, "chart_id": chart_id, "comparison": comp.star.auid, "check": check.star.auid, "observations": len(rows)}, indent=2))
    return output
