from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def bbox_to_wkt(min_lon: float, min_lat: float, max_lon: float, max_lat: float) -> str:
    return (
        "POLYGON(("
        f"{min_lon} {min_lat},"
        f"{max_lon} {min_lat},"
        f"{max_lon} {max_lat},"
        f"{min_lon} {max_lat},"
        f"{min_lon} {min_lat}"
        "))"
    )


def snap_gpt_path(explicit: str | None = None) -> str:
    if explicit:
        return explicit
    gpt = shutil.which("gpt")
    if not gpt:
        raise RuntimeError(
            "未找到 SNAP 的 `gpt` 命令。请先安装 ESA SNAP，并确保 `gpt` 在 PATH 中。"
        )
    return gpt


def preprocess_grd_zip_with_snap(
    input_zip: str | Path,
    output_tif: str | Path,
    *,
    bbox: tuple[float, float, float, float],
    graph_xml: str | Path,
    gpt: str | None = None,
) -> Path:
    input_zip = Path(input_zip)
    output_tif = Path(output_tif)
    output_tif.parent.mkdir(parents=True, exist_ok=True)
    graph_xml = Path(graph_xml)

    gpt = snap_gpt_path(gpt)
    wkt = bbox_to_wkt(*bbox)

    cmd = [
        gpt,
        str(graph_xml),
        f"-Pinput={input_zip}",
        f"-PsubsetWKT={wkt}",
        f"-Poutput={output_tif}",
    ]
    subprocess.run(cmd, check=True)
    return output_tif


def batch_preprocess(
    inputs: list[str | Path],
    outdir: str | Path,
    *,
    bbox: tuple[float, float, float, float],
    graph_xml: str | Path,
    gpt: str | None = None,
) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outs: list[Path] = []
    for inp in inputs:
        inp = Path(inp)
        out = outdir / (inp.stem + ".tif")
        outs.append(
            preprocess_grd_zip_with_snap(inp, out, bbox=bbox, graph_xml=graph_xml, gpt=gpt)
        )
    return outs

