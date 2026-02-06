"""
Sanchi oil spill (桑吉号) dataset builder.

Core idea:
- Use ASF Search API to reproducibly list Sentinel-1 GRD scenes intersecting an AOI and time window.
- Optionally download the GRD zips (Earthdata credentials typically required).
- Extract candidate oil-slick polygons from preprocessed GeoTIFFs (Sigma0 in dB recommended).
- Generate time-series metrics for numerical simulation validation.
"""

from .asf import search_scenes_geojson  # noqa: F401

