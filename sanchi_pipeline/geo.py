from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from pyproj import CRS, Geod, Transformer
from shapely.geometry import shape
from shapely.ops import transform, unary_union


WGS84 = CRS.from_epsg(4326)
GEOD = Geod(ellps="WGS84")


@dataclass(frozen=True)
class GeoMetrics:
    area_m2: float
    centroid_lon: float
    centroid_lat: float
    bbox_min_lon: float
    bbox_min_lat: float
    bbox_max_lon: float
    bbox_max_lat: float


def to_wgs84(geom, src_crs: CRS):
    if src_crs is None:
        return geom
    if CRS.from_user_input(src_crs) == WGS84:
        return geom
    tfm = Transformer.from_crs(src_crs, WGS84, always_xy=True)
    return transform(tfm.transform, geom)


def geodesic_area_m2(geom_wgs84) -> float:
    # pyproj expects GeoJSON-like mapping or shapely geometry; use geometry_area_perimeter.
    area, _ = GEOD.geometry_area_perimeter(geom_wgs84)
    return abs(float(area))


def metrics_from_wgs84(geom_wgs84) -> GeoMetrics:
    area_m2 = geodesic_area_m2(geom_wgs84)
    c = geom_wgs84.centroid
    minx, miny, maxx, maxy = geom_wgs84.bounds
    return GeoMetrics(
        area_m2=area_m2,
        centroid_lon=float(c.x),
        centroid_lat=float(c.y),
        bbox_min_lon=float(minx),
        bbox_min_lat=float(miny),
        bbox_max_lon=float(maxx),
        bbox_max_lat=float(maxy),
    )


def union_geoms(geoms: Iterable) -> object | None:
    geoms = [g for g in geoms if g and not g.is_empty]
    if not geoms:
        return None
    return unary_union(geoms)


def shape_from_geojson_geom(geojson_geom):
    return shape(geojson_geom)

