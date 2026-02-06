from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import requests

ASF_SEARCH_URL = "https://api.daac.asf.alaska.edu/services/search/param"


@dataclass(frozen=True)
class SearchParams:
    platforms: str  # comma-separated, e.g. "Sentinel-1A,Sentinel-1B"
    processing_level: str = "GRD_HD"
    start: str = "2018-01-06T00:00:00Z"
    end: str = "2018-01-20T23:59:59Z"
    bbox: tuple[float, float, float, float] = (123.5, 26.5, 126.5, 29.5)
    max_results: int = 2000


def search_scenes_geojson(params: SearchParams) -> dict[str, Any]:
    """
    Query ASF Search API and return a GeoJSON FeatureCollection.

    Notes:
    - This endpoint is publicly accessible for metadata queries.
    - For many AOI/time windows, results are small enough to fit in one request.
      We keep `max_results` configurable; pagination parameters are not reliably
      documented for all output formats.
    """
    min_lon, min_lat, max_lon, max_lat = params.bbox
    query = {
        "platform": params.platforms,
        "processingLevel": params.processing_level,
        "start": params.start,
        "end": params.end,
        "bbox": f"{min_lon},{min_lat},{max_lon},{max_lat}",
        "output": "geojson",
        "maxResults": str(params.max_results),
    }
    r = requests.get(ASF_SEARCH_URL, params=query, timeout=60)
    r.raise_for_status()
    return r.json()


def write_geojson(feature_collection: dict[str, Any], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(feature_collection, ensure_ascii=False, indent=2))
    return path


CSV_FIELDS = [
    "sceneName",
    "platform",
    "startTime",
    "stopTime",
    "flightDirection",
    "beamModeType",
    "polarization",
    "pathNumber",
    "frameNumber",
    "centerLon",
    "centerLat",
    "bytes",
    "processingLevel",
    "url",
    "browse",
    "md5sum",
]


def _get_prop(props: dict[str, Any], key: str) -> Any:
    v = props.get(key)
    if isinstance(v, list):
        return v[0] if v else ""
    return v if v is not None else ""


def features_to_rows(feature_collection: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for f in feature_collection.get("features", []):
        props = f.get("properties", {}) or {}
        row = {k: _get_prop(props, k) for k in CSV_FIELDS}
        rows.append(row)
    return rows


def write_csv(rows: Iterable[dict[str, Any]], path: str | Path, fields: list[str] | None = None) -> Path:
    import csv

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or CSV_FIELDS
    with path.open("w", newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})
    return path

