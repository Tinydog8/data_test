from __future__ import annotations

import argparse
import glob
from pathlib import Path

from .asf import SearchParams, features_to_rows, search_scenes_geojson, write_csv, write_geojson
from .defaults import (
    DEFAULT_BBOX,
    DEFAULT_END,
    DEFAULT_PLATFORMS,
    DEFAULT_PROCESSING_LEVEL,
    DEFAULT_START,
)
from .download import download_from_geojson
from .extract import batch_extract
from .snap import batch_preprocess
from .timeseries import polygons_dir_to_timeseries_csv


def _bbox_from_args(vals: list[str]) -> tuple[float, float, float, float]:
    if len(vals) != 4:
        raise SystemExit("`--bbox` 需要 4 个数字：min_lon min_lat max_lon max_lat")
    return tuple(map(float, vals))  # type: ignore[return-value]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="sanchi_pipeline", description="桑吉号溢油数据集构建（可复现）")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_search = sub.add_parser("search", help="检索 Sentinel-1 GRD 场景清单（ASF Search）")
    p_search.add_argument("--platforms", default=DEFAULT_PLATFORMS, help="例如 Sentinel-1A,Sentinel-1B")
    p_search.add_argument("--processing-level", default=DEFAULT_PROCESSING_LEVEL)
    p_search.add_argument("--start", default=DEFAULT_START)
    p_search.add_argument("--end", default=DEFAULT_END)
    p_search.add_argument("--bbox", nargs=4, default=list(DEFAULT_BBOX))
    p_search.add_argument("--max-results", type=int, default=2000)
    p_search.add_argument("--outdir", default="data/scenes")

    p_dl = sub.add_parser("download", help="按场景清单下载 GRD ZIP（可能需要 Earthdata）")
    p_dl.add_argument("--scenes", required=True, help="`search` 输出的 scenes.geojson")
    p_dl.add_argument("--outdir", default="data/raw_zips")
    p_dl.add_argument("--username", default=None)
    p_dl.add_argument("--password", default=None)
    p_dl.add_argument("--overwrite", action="store_true")
    p_dl.add_argument("--limit", type=int, default=None)

    p_snap = sub.add_parser("snap-preprocess", help="用 SNAP gpt 预处理 GRD ZIP → GeoTIFF")
    p_snap.add_argument("--input-zip", nargs="+", required=True, help="一个或多个 .zip（支持通配符）")
    p_snap.add_argument("--outdir", default="data/geotiff")
    p_snap.add_argument("--bbox", nargs=4, default=list(DEFAULT_BBOX))
    p_snap.add_argument("--graph-xml", default="snap_graphs/s1_grd_preprocess_subset_tc.xml")
    p_snap.add_argument("--gpt", default=None, help="gpt 可执行文件路径（可选）")

    p_ex = sub.add_parser("extract", help="从预处理 GeoTIFF 提取油膜候选多边形（GeoJSON）")
    p_ex.add_argument("--input", required=True, help="GeoTIFF 文件或目录")
    p_ex.add_argument("--outdir", default="data/oil_polygons")
    p_ex.add_argument("--band", type=int, default=1)
    p_ex.add_argument("--input-is-db", action="store_true", default=True)
    p_ex.add_argument("--input-is-linear", action="store_true", default=False, help="输入是线性 sigma0（会转 dB）")
    p_ex.add_argument("--threshold-k", type=float, default=2.5)
    p_ex.add_argument("--threshold-quantile", type=float, default=0.02, help="覆盖率过大时的回退分位数阈值")
    p_ex.add_argument("--max-coverage", type=float, default=0.2, help="阈值分割后 True 覆盖率上限，超过则回退到分位数阈值")
    p_ex.add_argument("--min-area-m2", type=float, default=50_000.0)
    p_ex.add_argument("--open-px", type=int, default=2)
    p_ex.add_argument("--close-px", type=int, default=4)
    p_ex.add_argument("--fill-holes", action="store_true", help="是否填充油膜掩膜孔洞（默认关闭以避免误填充）")

    p_ts = sub.add_parser("timeseries", help="从 GeoJSON 多边形输出时间序列 CSV")
    p_ts.add_argument("--polygons", required=True, help="extract 输出目录")
    p_ts.add_argument("--out", required=True, help="输出 CSV 路径")

    args = p.parse_args(argv)

    if args.cmd == "search":
        bbox = _bbox_from_args(args.bbox)
        params = SearchParams(
            platforms=args.platforms,
            processing_level=args.processing_level,
            start=args.start,
            end=args.end,
            bbox=bbox,
            max_results=args.max_results,
        )
        fc = search_scenes_geojson(params)
        outdir = Path(args.outdir)
        gj = write_geojson(fc, outdir / "scenes.geojson")
        rows = features_to_rows(fc)
        csv_path = write_csv(rows, outdir / "scenes.csv")
        print(f"已输出: {gj}")
        print(f"已输出: {csv_path}")
        print(f"场景数: {len(rows)}")
        return 0

    if args.cmd == "download":
        written = download_from_geojson(
            args.scenes,
            args.outdir,
            username=args.username,
            password=args.password,
            overwrite=args.overwrite,
            limit=args.limit,
        )
        print(f"下载完成: {len(written)} 个文件 → {Path(args.outdir)}")
        return 0

    if args.cmd == "snap-preprocess":
        bbox = _bbox_from_args(args.bbox)
        inputs: list[str] = []
        for pat in args.input_zip:
            matches = glob.glob(pat)
            if matches:
                inputs.extend(matches)
            else:
                inputs.append(pat)
        outs = batch_preprocess(
            inputs,
            args.outdir,
            bbox=bbox,
            graph_xml=args.graph_xml,
            gpt=args.gpt,
        )
        print(f"预处理完成: {len(outs)} 个 GeoTIFF → {Path(args.outdir)}")
        return 0

    if args.cmd == "extract":
        input_is_db = True
        if args.input_is_linear:
            input_is_db = False
        outs = batch_extract(
            args.input,
            args.outdir,
            band=args.band,
            input_is_db=input_is_db,
            threshold_k=args.threshold_k,
            threshold_quantile=args.threshold_quantile,
            max_coverage=args.max_coverage,
            min_area_m2=args.min_area_m2,
            open_px=args.open_px,
            close_px=args.close_px,
            fill_holes=args.fill_holes,
        )
        print(f"提取完成: {len(outs)} 个 GeoJSON → {Path(args.outdir)}")
        return 0

    if args.cmd == "timeseries":
        out = polygons_dir_to_timeseries_csv(args.polygons, args.out)
        print(f"已输出: {out}")
        return 0

    raise SystemExit("未知命令")


if __name__ == "__main__":
    raise SystemExit(main())

