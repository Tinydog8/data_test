from __future__ import annotations

import re
from datetime import datetime, timezone


_S1_NAME_RE = re.compile(
    r"^(?P<platform>S1[AB])_.*?_(?P<start>\d{8}T\d{6})_(?P<stop>\d{8}T\d{6})_"
)


def parse_s1_times_from_name(scene_name: str) -> tuple[str | None, str | None]:
    """
    Parse start/stop timestamps from a Sentinel-1 scene name.

    Example:
    S1A_IW_GRDH_1SDV_20180117T095411_20180117T095436_020193_022732_B76C
    """
    m = _S1_NAME_RE.match(scene_name)
    if not m:
        return None, None

    def _to_iso(s: str) -> str:
        dt = datetime.strptime(s, "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
        return dt.isoformat().replace("+00:00", "Z")

    return _to_iso(m.group("start")), _to_iso(m.group("stop"))

