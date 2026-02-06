from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .geo import metrics_from_wgs84, shape_from_geojson_geom, union_geoms
from .utils import parse_s1_times_from_name


def _scene_from_featurecollection(fc: dict[str, Any], fallback_name: str) -> tuple[str, str | None, str | None]:
    props = (fc.get("properties") or {}) if isinstance(fc, dict) else {}
    scene = props.get("sceneName") or fc.get("name") or fallback_name
    start = props.get("startTime")
    stop = props.get("stopTime")
    if not start or not stop:
        s, e = parse_s1_times_from_name(scene)
        start = start or s
        stop = stop or e
    return str(scene), start, stop


def polygons_dir_to_timeseries_csv(polygons_dir: str | Path, out_csv: str | Path) -> Path:
    polygons_dir = Path(polygons_dir)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for p in sorted(polygons_dir.rglob("*.geojson")):
        fc = json.loads(p.read_text(encoding="utf-8"))
        scene, start, stop = _scene_from_featurecollection(fc, fallback_name=p.stem)

        geoms = []
        for feat in fc.get("features", []):
            g = shape_from_geojson_geom(feat.get("geometry"))
            if g and not g.is_empty:
                geoms.append(g)
        u = union_geoms(geoms)
        if u is None:
            continue
        m = metrics_from_wgs84(u)
        rows.append(
            {
                "sceneName": scene,
                "startTime": start or "",
                "stopTime": stop or "",
                "area_m2": m.area_m2,
                "centroid_lon": m.centroid_lon,
                "centroid_lat": m.centroid_lat,
                "bbox_min_lon": m.bbox_min_lon,
                "bbox_min_lat": m.bbox_min_lat,
                "bbox_max_lon": m.bbox_max_lon,
                "bbox_max_lat": m.bbox_max_lat,
                "source_geojson": str(p),
            }
        )

    rows.sort(key=lambda r: (r.get("startTime", ""), r.get("sceneName", "")))

    fields = [
        "sceneName",
        "startTime",
        "stopTime",
        "area_m2",
        "centroid_lon",
        "centroid_lat",
        "bbox_min_lon",
        "bbox_min_lat",
        "bbox_max_lon",
        "bbox_max_lat",
        "source_geojson",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    return out_csv

