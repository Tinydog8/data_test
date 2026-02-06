from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Iterable

import requests
from tqdm import tqdm


def iter_scene_downloads(feature_collection: dict[str, Any]) -> Iterable[tuple[str, str, int | None]]:
    for f in feature_collection.get("features", []):
        props = f.get("properties", {}) or {}
        url = props.get("url") or ""
        name = props.get("fileName") or props.get("sceneName") or ""
        size = props.get("bytes")
        try:
            size_int = int(size) if size is not None else None
        except Exception:
            size_int = None
        if url and name:
            yield url, name, size_int


def _env_or_none(key: str) -> str | None:
    v = os.environ.get(key)
    return v if v else None


def download_from_geojson(
    scenes_geojson: str | Path,
    outdir: str | Path,
    username: str | None = None,
    password: str | None = None,
    overwrite: bool = False,
    limit: int | None = None,
) -> list[Path]:
    import json

    scenes_geojson = Path(scenes_geojson)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fc = json.loads(scenes_geojson.read_text(encoding="utf-8"))

    username = username or _env_or_none("EARTHDATA_USERNAME")
    password = password or _env_or_none("EARTHDATA_PASSWORD")

    session = requests.Session()
    if username and password:
        session.auth = (username, password)

    downloads = list(iter_scene_downloads(fc))
    if limit is not None:
        downloads = downloads[: max(0, int(limit))]

    written: list[Path] = []
    for url, name, expected_size in tqdm(downloads, desc="download", unit="file"):
        out = outdir / name
        if out.exists() and not overwrite:
            if expected_size is None or out.stat().st_size == expected_size:
                written.append(out)
                continue

        with session.get(url, stream=True, timeout=120, allow_redirects=True) as r:
            r.raise_for_status()
            total = int(r.headers.get("Content-Length", "0") or "0") or expected_size or 0
            tmp = out.with_suffix(out.suffix + ".part")
            with tmp.open("wb") as fp, tqdm(
                total=total if total > 0 else None,
                desc=out.name,
                unit="B",
                unit_scale=True,
                leave=False,
            ) as bar:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if not chunk:
                        continue
                    fp.write(chunk)
                    bar.update(len(chunk))
            tmp.replace(out)
        written.append(out)

    return written

