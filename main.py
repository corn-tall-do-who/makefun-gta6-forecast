"""GTA VI Twitch CCU forecaster, variant B (corn-tall-do-who).

Log-scale quantile forecaster anchored on the first hour's print.

Model, in prose. The first hour's viewer count is the anchor: every later
hour is forecast as a multiple of it. The multiple's centre (the median) is
a piecewise-linear function of the hour in log space, with three segments
-- the early hours, the middle of the day and the late hours -- each with
its own slope. The slopes are fitted at run time from the past big launches
the case itself carries: for every comparable launch, take its release-day
hourly curve, divide each hour by that launch's hour-1 count, and take the
log. Only launches whose hour-1 count is within a factor of five of this
case's hour-1 count are used as peers -- a launch that had barely started
at 00:00 UTC behaves nothing like one that was already at full launch
level, and the anchor's size is what tells them apart. When too few peers
qualify, the band widens by a factor of five at a time (25x, then 125x)
before falling back to every comparable with a usable release day.

The spread widens with the hour. Around the fitted median, the residuals
of the peer launches are collected over a five-hour window centred on the
hour being forecast, and each level's offset is that window's empirical
quantile of the residuals -- so the 5% and 95% bands are honest about how
far peer launches actually strayed, including the long upper tail of
launches that unlock later in the day, instead of assuming a symmetric
bell shape.

Everything is fitted from the case file at run time, offline, with the
standard library only. No randomness anywhere: the same input bytes give
the same output bytes. If the comparables cannot support a fit at all, the
forecast falls back to holding the hour-1 print with a spread that widens
like a random walk, so the output is always a valid forecast.
"""

import json
import math
import os
from datetime import datetime
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

INPUT_DIR = Path(os.environ.get("MAKEFUN_INPUT_DIR", "/input"))
OUTPUT_DIR = Path(os.environ.get("MAKEFUN_OUTPUT_DIR", "/output"))

# The quantile levels the contest scores, in the order each row lists them.
LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
# Fallback spread: log-units per square-root hour, used only when no fit
# is possible at all (mirrors the house baseline's widening, flat median).
FALLBACK_SIGMA = 0.15
FALLBACK_Z = [-1.645, -0.674, 0.0, 0.674, 1.645]
# Peer band: a comparable counts as a peer when its hour-1 count is within
# this factor of the case's hour-1 count, either way.
PEER_FACTOR = 5.0
MIN_PEERS = 5
# A comparable is usable when its release day holds at least this many
# of the 24 hourly points including hour 1.
MIN_DAY_HOURS = 20
# Half-width, in hours, of the window pooling residuals for one hour.
SPREAD_WINDOW = 2
# Breakpoints of the three median segments, in x = hour - 1.
BREAK_1 = 7
BREAK_2 = 15


def parse_ts(text):
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def release_day(comp):
    """{hour: viewers} for hours 1..24 of a comparable's release day."""
    start = parse_ts(comp["day_start"])
    out = {}
    for point in comp.get("hourly") or []:
        viewers = point.get("viewers")
        if not viewers or viewers <= 0:
            continue
        delta = (parse_ts(point["from"]) - start).total_seconds()
        if 0 <= delta < 24 * 3600 and delta % 3600 == 0:
            out[int(delta // 3600) + 1] = viewers
    return out


def quantile(sorted_values, tau):
    """Linear-interpolated tau-quantile of an already sorted list."""
    n = len(sorted_values)
    if n == 0:
        return 0.0
    if n == 1:
        return sorted_values[0]
    pos = tau * (n - 1)
    low = int(math.floor(pos))
    high = min(low + 1, n - 1)
    frac = pos - low
    return sorted_values[low] * (1.0 - frac) + sorted_values[high] * frac


def solve_linear(matrix, rhs):
    """Small dense solve by Gauss elimination with partial pivoting."""
    n = len(rhs)
    rows = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        piv = rows[col][col]
        if abs(piv) < 1e-12:
            continue
        for r in range(n):
            if r != col and rows[r][col] != 0.0:
                factor = rows[r][col] / piv
                for c in range(col, n + 1):
                    rows[r][c] -= factor * rows[col][c]
    out = []
    for i in range(n):
        diag = rows[i][i]
        out.append(rows[i][n] / diag if abs(diag) > 1e-12 else 0.0)
    return out


def least_squares(design, targets):
    """Least-squares coefficients for a small design matrix."""
    k = len(design[0])
    ata = [[0.0] * k for _ in range(k)]
    atb = [0.0] * k
    for row, y in zip(design, targets):
        for i in range(k):
            atb[i] += row[i] * y
            for j in range(k):
                ata[i][j] += row[i] * row[j]
    return solve_linear(ata, atb)


def segment_basis(x):
    """Continuous piecewise-linear basis in x = hour - 1.

    [1, early part, middle part, late part]: the first slope runs to hour 8,
    the second from hour 8 to hour 16, the third after hour 16.
    """
    return [
        1.0,
        min(x, BREAK_1),
        min(max(x - BREAK_1, 0), BREAK_2 - BREAK_1),
        max(x - BREAK_2, 0),
    ]


def median_curve(peers):
    """Per-hour median log ratio to hour 1, over the peer launches."""
    medians = {}
    for hour in range(2, 25):
        ratios = sorted(
            math.log(day[hour] / day[1]) for day in peers if hour in day
        )
        if ratios:
            medians[hour] = quantile(ratios, 0.5)
    return medians


def fit_median(peers):
    """Piecewise-linear median log ratio, fitted to the per-hour medians.

    Returns the four coefficients, or None when too few medians exist.
    """
    medians = median_curve(peers)
    hours = sorted(medians)
    if len(hours) < 4:
        return None
    coef = least_squares(
        [segment_basis(h - 1) for h in hours], [medians[h] for h in hours]
    )
    return coef


def median_at(coef, hour):
    total = 0.0
    for a, b in zip(coef, segment_basis(hour - 1)):
        total += a * b
    return total


def residual_quantiles(peers, coef, levels):
    """Empirical residual quantile offsets per hour, per level.

    For hour h, pool the peers' residuals (log ratio minus fitted median)
    over hours h-2..h+2 and take each level's quantile of that pool.
    """
    offsets = {}
    for hour in range(2, 25):
        lo = max(2, hour - SPREAD_WINDOW)
        hi = min(24, hour + SPREAD_WINDOW)
        pooled = []
        for near in range(lo, hi + 1):
            centre = median_at(coef, near)
            for day in peers:
                if near in day:
                    pooled.append(math.log(day[near] / day[1]) - centre)
        if len(pooled) < MIN_PEERS:
            pooled = []
            for near in range(2, 25):
                centre = median_at(coef, near)
                for day in peers:
                    if near in day:
                        pooled.append(math.log(day[near] / day[1]) - centre)
        pooled.sort()
        offsets[hour] = [quantile(pooled, tau) for tau in levels]
    return offsets


def round2(value):
    """Round to 2 decimals, half up; whole numbers come back as ints."""
    rounded = Decimal(repr(value)).quantize(
        Decimal("0.01"), rounding=ROUND_HALF_UP
    )
    if rounded == rounded.to_integral_value():
        return int(rounded)
    return float(rounded)


def clamp_row(values):
    """Enforce finite, positive, capped, non-decreasing quantiles."""
    fixed = []
    prev = 1.0
    for value in values:
        if not math.isfinite(value):
            value = prev
        value = min(max(value, 1.0), 1e9)
        value = max(value, prev)
        fixed.append(round2(value))
        prev = fixed[-1]
    return fixed


def flat_forecast(anchor, hours, levels):
    """Fallback: hold the anchor, spread widening like a random walk."""
    rows = []
    for hour in hours:
        spread = FALLBACK_SIGMA * math.sqrt(hour - 1)
        qs = [
            anchor * math.exp(z * spread)
            for z in FALLBACK_Z[: len(levels)]
        ]
        rows.append({"hour": hour, "q": clamp_row(qs)})
    return {"hours": rows}


def forecast_case(case):
    levels = list(case.get("levels") or LEVELS)
    hours = [row["hour"] for row in case["forecast"]]
    first = next(
        (k for k in case.get("known") or [] if k.get("hour") == 1), None
    )
    anchor = 1.0
    if first is not None and first.get("viewers"):
        anchor = max(1.0, float(first["viewers"]))
    try:
        days = []
        for comp in sorted(
            case.get("comparables") or [], key=lambda c: str(c.get("name"))
        ):
            day = release_day(comp)
            if 1 in day and len(day) >= MIN_DAY_HOURS:
                days.append(day)
        if not days:
            raise ValueError("no usable comparables")
        peers = []
        for factor in (PEER_FACTOR, PEER_FACTOR**2, PEER_FACTOR**3):
            peers = [
                d
                for d in days
                if anchor / factor <= d[1] <= anchor * factor
            ]
            if len(peers) >= MIN_PEERS:
                break
        if len(peers) < MIN_PEERS:
            peers = days
        coef = fit_median(peers)
        if coef is None:
            raise ValueError("median fit failed")
        offsets = residual_quantiles(peers, coef, levels)
        rows = []
        for hour in hours:
            centre = median_at(coef, hour)
            qs = [
                anchor * math.exp(centre + off)
                for off in offsets[hour]
            ]
            rows.append({"hour": hour, "q": clamp_row(qs)})
        return {"hours": rows}
    except Exception:
        return flat_forecast(anchor, hours, levels)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(INPUT_DIR.glob("*.in")):
        case = json.loads(path.read_text(encoding="utf-8"))
        forecast = forecast_case(case)
        out_path = OUTPUT_DIR / (path.stem + ".out")
        out_path.write_text(
            json.dumps(forecast, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
