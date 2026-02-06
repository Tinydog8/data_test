from __future__ import annotations

import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from shapely.geometry import mapping, shape

from .geo import metrics_from_wgs84, to_wgs84, union_geoms
from .utils import parse_s1_times_from_name


def _pixel_area_m2(transform, crs, mid_lat: float) -> float:
    # If projected (meters), use affine determinant.
    try:
        from pyproj import CRS

        if crs and not CRS.from_user_input(crs).is_geographic:
            return abs(float(transform.a * transform.e))
    except Exception:
        pass

    # Geographic degrees: approximate.
    # transform.a ~ lon degrees/pixel, transform.e ~ lat degrees/pixel (negative).
    res_lon = abs(float(transform.a))
    res_lat = abs(float(transform.e))
    m_per_deg_lat = 110_574.0
    m_per_deg_lon = 111_320.0 * math.cos(math.radians(mid_lat))
    return (m_per_deg_lon * res_lon) * (m_per_deg_lat * res_lat)


def _robust_threshold_db(img_db: np.ndarray, k: float = 2.5) -> float:
    v = img_db[np.isfinite(img_db)]
    if v.size == 0:
        return float("nan")
    med = float(np.nanmedian(v))
    mad = float(np.nanmedian(np.abs(v - med)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 1e-6:
        sigma = float(np.nanstd(v)) or 1.0
    return med - k * sigma


def _quantile_threshold_db(img_db: np.ndarray, q: float = 0.02) -> float:
    v = img_db[np.isfinite(img_db)]
    if v.size == 0:
        return float("nan")
    q = float(q)
    q = min(max(q, 0.0), 1.0)
    return float(np.nanquantile(v, q))


def extract_oil_polygons_geotiff(
    geotiff: str | Path,
    out_geojson: str | Path,
    *,
    band: int = 1,
    input_is_db: bool = True,
    threshold_k: float = 2.5,
    threshold_quantile: float = 0.02,
    max_coverage: float = 0.2,
    min_area_m2: float = 50_000.0,
    open_px: int = 2,
    close_px: int = 4,
    fill_holes: bool = False,
) -> Path:
    """
    Extract candidate oil-slick polygons from a preprocessed GeoTIFF.

    Assumptions:
    - GeoTIFF is already calibrated / terrain-corrected.
    - Pixel values are Sigma0 in dB by default (`input_is_db=True`).
    """
    try:
        import rasterio
        import rasterio.features
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "缺少依赖：需要安装 rasterio 才能从 GeoTIFF 提取多边形。请先 `pip install -r requirements.txt`。"
        ) from e

    try:
        from scipy import ndimage as ndi
    except Exception as e:  # pragma: no cover
        raise RuntimeError("缺少依赖：需要安装 scipy。") from e

    geotiff = Path(geotiff)
    out_geojson = Path(out_geojson)
    out_geojson.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(geotiff) as ds:
        arr = ds.read(band).astype("float32")
        nodata = ds.nodata
        if nodata is not None:
            arr[arr == nodata] = np.nan
        if ds.mask_flag_enums and any(ds.mask_flag_enums):
            # Apply dataset mask if present.
            mask = ds.dataset_mask() == 0
            arr[mask] = np.nan

        if not input_is_db:
            arr = 10.0 * np.log10(np.maximum(arr, 1e-10))

        valid = np.isfinite(arr)
        thr = _robust_threshold_db(arr, k=threshold_k)
        oil = valid & (arr < thr)
        coverage = float(oil.sum()) / float(max(1, valid.sum()))
        if coverage > max_coverage:
            # Robust threshold may be too loose for some scenes; fall back to a fixed dark-tail quantile.
            thr = _quantile_threshold_db(arr, q=threshold_quantile)
            oil = valid & (arr < thr)

        if open_px > 0:
            oil = ndi.binary_opening(oil, iterations=int(open_px))
        if close_px > 0:
            oil = ndi.binary_closing(oil, iterations=int(close_px))
        if fill_holes:
            oil = ndi.binary_fill_holes(oil)

        # Remove small components
        labeled, n = ndi.label(oil)
        if n > 0:
            # approx pixel area
            mid_lat = float(ds.bounds.bottom + ds.bounds.top) / 2.0
            pix_area = _pixel_area_m2(ds.transform, ds.crs, mid_lat=mid_lat)
            min_pixels = int(max(1, math.ceil(min_area_m2 / max(pix_area, 1e-6))))
            counts = np.bincount(labeled.ravel())
            keep = counts >= min_pixels
            # label=0 is background, never keep
            if keep.size > 0:
                keep[0] = False
            oil = keep[labeled]

        # Important: mask should represent "valid pixels", not the oil mask itself.
        # If we mask out the background (0) pixels, `shapes()` may effectively dissolve
        # disconnected oil pixels into one big polygon envelope.
        shapes = list(
            rasterio.features.shapes(oil.astype("uint8"), mask=valid, transform=ds.transform)
        )
        geoms = []
        for geom_mapping, val in shapes:
            if int(val) != 1:
                continue
            g = to_wgs84(
                shape(geom_mapping),
                ds.crs,
            )
            if g.is_empty:
                continue
            geoms.append(g)

    union = union_geoms(geoms)
    scene_name = geotiff.stem
    start_time, stop_time = parse_s1_times_from_name(scene_name)

    features: list[dict[str, Any]] = []
    total_area = 0.0
    for i, g in enumerate(geoms):
        m = metrics_from_wgs84(g)
        total_area += m.area_m2
        features.append(
            {
                "type": "Feature",
                "id": f"{scene_name}:{i}",
                "geometry": mapping(g),
                "properties": {
                    "sceneName": scene_name,
                    "startTime": start_time,
                    "stopTime": stop_time,
                    "area_m2": m.area_m2,
                    "centroid_lon": m.centroid_lon,
                    "centroid_lat": m.centroid_lat,
                    "bbox_min_lon": m.bbox_min_lon,
                    "bbox_min_lat": m.bbox_min_lat,
                    "bbox_max_lon": m.bbox_max_lon,
                    "bbox_max_lat": m.bbox_max_lat,
                    "source_raster": str(geotiff),
                },
            }
        )

    collection = {
        "type": "FeatureCollection",
        "name": scene_name,
        "features": features,
        "properties": {
            "sceneName": scene_name,
            "startTime": start_time,
            "stopTime": stop_time,
            "threshold_k": threshold_k,
            "threshold_quantile": threshold_quantile,
            "max_coverage": max_coverage,
            "min_area_m2": min_area_m2,
            "open_px": open_px,
            "close_px": close_px,
            "fill_holes": fill_holes,
            "total_area_m2": float(total_area),
            "union_area_m2": metrics_from_wgs84(union).area_m2 if union else 0.0,
            "union_metrics": asdict(metrics_from_wgs84(union)) if union else None,
        },
    }

    out_geojson.write_text(json.dumps(collection, ensure_ascii=False, indent=2), encoding="utf-8")
    return out_geojson


def batch_extract(
    input_path: str | Path,
    outdir: str | Path,
    *,
    glob: str = "*.tif",
    **kwargs,
) -> list[Path]:
    input_path = Path(input_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if input_path.is_dir():
        patterns = [glob]
        # common alias: *.tiff
        if glob.endswith("*.tif"):
            patterns.append("*.tiff")
        files: list[Path] = []
        for pat in patterns:
            files.extend(input_path.rglob(pat))
        files = sorted(set(files))
    else:
        files = [input_path]

    outputs: list[Path] = []
    for f in files:
        out = outdir / (f.stem + ".geojson")
        outputs.append(extract_oil_polygons_geotiff(f, out, **kwargs))
    return outputs

