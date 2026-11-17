#!/usr/bin/env python3
"""Forecast GTA VI's launch-day Twitch viewers, hour by hour (hours 2-24).

make.fun competition gta6-twitch-ccu-forecast-01, round v30. Python 3.13,
standard library only, no network. Reads one case per *.in file in
$MAKEFUN_INPUT_DIR (else /input/) and writes one forecast per case to
$MAKEFUN_OUTPUT_DIR (else /output/): {"hours": [{"hour": h,
"q": [q05, q25, q50, q75, q95]}, ...]}. Score is pinball loss on the log
scale; the program outputs symmetric-lognormal quantiles around a
log-median with a per-side log spread.

Model (v30 packaged cell; v29 values in parentheses where changed)
------------------------------------------------
Four member predictors, each a log-space blend of one per-hour OLS view
and two one-factor shape views (console-only and all launches):

* Y5   -- shape pools: tier != "ea"; shape regressions on [lnV1, g12];
  population sd; mix 0.5/0.2/0.3; spread = sqrt of the weighted mean
  of the squared view spreads, floored at sqrt(0.02).
* Z0   -- same as Y5, but the naive early-access filter ("early access"
  as a substring of the unlock text) instead of the tier filter.
* Y1   -- same pools/spreads as Z0, but its shape (level, amplitude)
  regressions run on [lnV1, lnPREM] -- the v1 original features --
  instead of [lnV1, g12].
* clone -- the v2 reconstruction: hyphen early-access filter
  ("early-access" substring); shape regressions on [lnV1, lnPREM, g12];
  sample sd; whole-library shape pools (no +-2.5 window); mix
  0.6/0.2/0.2; spread = weighted mean of view spreads x 0.81.

The OLS view is shared by all members: per hour, ln(viewers) on
[lnV1, lnPREM, g12] over launches within +-2.5 log units of the
target's (clamped) lnV1, whole library when the window holds fewer
than 10 launches, median-recentered, spread = residual sd x 0.8.

A shape view: the median log-ratio curve m_h = per-hour median of
ln(V_h / V_1) over its pool; per launch a level/amplitude fit
r_h ~= c + s * m_h (launches with >= 20 hours); c and s regressed on
the member's features across launches (>= 5 fits); amplitude clamped
to [0, 2]; the view's median is lnV1 + c* + s* * m_h and its spread
is the per-hour sd of the fit residuals x 1.0. The target's size fed
to every regression is clamped to the library's lnV1 range plus 0.5
log units of headroom; the level anchor stays the true hour-1 print.

Combination (per hour h, bands early h2-6 / mid h7-14 / late h15-24):
median source early = Y1, mid = 0.5*Y5 + 0.5*clone, late =
median(Y5, Z0, Y1); ramps early 0.00 -> +0.03 (was +0.025), mid
0.00 -> -0.12 (was -0.115), late 0.00 -> +0.05; size shift
-0.02 * max(0, lnV1 - 11.5) (was -0.015); disagreement shift
(0.40 + 0.10 * max(0, lnV1 - 11.5)) * sd of the four member medians
(was 0.35 + 0.05 * size), except the late band, which uses 0.35 as
its base coefficient; spread = 0.975 * (0.45*sig(Y5) + 0.35*sig(Z0)
+ 0.20*sig(clone)), times (1 + 0.80 * sd) (was 0.75), with the 0.975
replaced in late hours by a ramp 0.975 -> 0.90 (was 0.92);
per-band asymmetric scales early (0.90, 1.10), mid (1.10, 0.80),
late (0.95, 0.85); both side spreads capped at 0.315.

Honest caveat: every band parameter was chosen on the four practice
launches and on leave-one-out over the comparables across many rounds
of selection, so neither gate is a clean holdout and the margins are
thin and selection-biased. The disagreement shift is an empirical
regularity, not a mechanism, and the final launch (lnV1 ~ 14+) sits
beyond the library's largest launch, where the shift extrapolates.
Treat this as the best verified fallback, not a guarantee.
"""

import json
import math
import os
import sys
from datetime import datetime, timedelta, timezone

Z05 = 1.6448536269514722
Z25 = 0.6744897501960817
HOURS = list(range(2, 25))
EARLY = list(range(2, 7))
MID = list(range(7, 15))
LATE = list(range(15, 25))
EPS = 1e-9

# Packaged-cell parameters (v29 baseline). The v30 round tunes these.
P = {
    "early_ramp": 0.03,
    "mid_ramp": -0.12,
    "late_ramp": 0.05,
    "size_slope": -0.02,
    "size_pivot": 11.5,
    "shift_base": 0.40,
    "shift_slope": 0.10,
    "spread_base": 0.975,
    "spread_late_end": 0.90,
    "disag_k": 0.8,
    "spread_cap": 0.315,
    "asym": {"early": (0.90, 1.10), "mid": (1.10, 0.80), "late": (0.95, 0.85)},
    "spread_w": (0.45, 0.35, 0.20),
    # v30 sweep knobs (defaults reproduce the v29 cell exactly).
    "shift_late": 0.35,          # late-band shift coefficient override
    "early_src": "y1",           # "y1" | "med3" | "y1z0"
    "mid_src": "y5clone",        # "y5clone" | "med4" | "y5clonez0"
    "late_src": "med3",          # "med3" | "med4" | "mean3"
    "spread_size_slope": 0.0,    # spread *= 1 + slope * max(0, lnV1-11.5)
}


# ---------------------------------------------------------------------------
# Small numeric helpers
# ---------------------------------------------------------------------------

def _med(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return float("nan")
    if n % 2 == 1:
        return s[n // 2]
    return 0.5 * (s[n // 2 - 1] + s[n // 2])


def _pstdev(xs):
    if len(xs) < 2:
        return 0.0
    m = sum(xs) / len(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / len(xs))


def _stdev(xs):
    if len(xs) < 2:
        return 0.0
    m = sum(xs) / len(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def _solve(mat, vec):
    """Gaussian elimination with partial pivoting; None when singular."""
    n = len(vec)
    a = [row[:] + [vec[i]] for i, row in enumerate(mat)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-12:
            return None
        a[col], a[piv] = a[piv], a[col]
        pv = a[col][col]
        for r in range(n):
            if r != col and a[r][col] != 0.0:
                f = a[r][col] / pv
                for c in range(col, n + 1):
                    a[r][c] -= f * a[col][c]
    return [a[i][n] / a[i][i] for i in range(n)]


def _lstsq(rows, ys):
    """Least squares for rows (lists, intercept included). None if degenerate."""
    if not rows:
        return None
    k = len(rows[0])
    ata = [[0.0] * k for _ in range(k)]
    aty = [0.0] * k
    for row, y in zip(rows, ys):
        for i in range(k):
            aty[i] += row[i] * y
            for j in range(k):
                ata[i][j] += row[i] * row[j]
    for i in range(k):
        ata[i][i] += 1e-10
    return _solve(ata, aty)


def _parse_ts(text):
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def _band(hour):
    if hour <= 6:
        return "early"
    if hour <= 14:
        return "mid"
    return "late"


def _ramp_at(hour, start, end, first, last):
    if last <= first:
        return end
    frac = (hour - first) / (last - first)
    return start + (end - start) * frac


# ---------------------------------------------------------------------------
# Case loading and feature extraction
# ---------------------------------------------------------------------------

def _hourly_map(points):
    out = {}
    for p in points or []:
        try:
            out[_parse_ts(p["from"])] = float(p["viewers"])
        except (KeyError, TypeError, ValueError):
            continue
    return out


def _tier_of(unlock):
    """Tier from the 'Tier: X' annotation in the unlock text."""
    text = unlock or ""
    low = text.lower()
    marker = "tier:"
    idx = low.find(marker)
    if idx < 0:
        return "comparable"
    rest = low[idx + len(marker):].strip()
    token = rest.split()[0].strip(".,;") if rest.split() else "comparable"
    if token.startswith("early-access"):
        return "ea"
    if token == "weak":
        return "weak"
    if token.startswith("target-history"):
        return "target-history"
    return "comparable"


def _launch_features(comp):
    """Features of one comparable launch, or None when unusable."""
    try:
        day = _parse_ts(comp["day_start"])
    except (KeyError, TypeError, ValueError):
        return None
    hourly = _hourly_map(comp.get("hourly") or [])
    v1 = hourly.get(day)
    if not v1 or v1 <= 0:
        return None
    ln_v1 = math.log(v1)
    pre = [v for t, v in hourly.items()
           if day - timedelta(days=7) <= t < day and v > 0]
    if not pre:
        return None
    ln_prem = math.log(sum(pre) / len(pre))
    point12 = hourly.get(day - timedelta(hours=12))
    g12 = None
    if point12 and point12 > 0:
        g12 = ln_v1 - math.log(point12)
    ratios = {}
    for h in HOURS:
        v = hourly.get(day + timedelta(hours=h - 1))
        if v and v > 0:
            ratios[h] = math.log(v / v1)
    unlock = comp.get("unlock") or ""
    return {
        "name": comp.get("name", ""),
        "ln_v1": ln_v1,
        "ln_prem": ln_prem,
        "g12": g12,
        "console": "PC" not in (comp.get("platforms") or ""),
        "tier": _tier_of(unlock),
        "unlock_low": unlock.lower(),
        "ratios": ratios,
        "hourly": hourly,
        "day": day,
    }


def _target_features(case):
    """Features of the case's own launch (the one being forecast)."""
    day = _parse_ts(case["day_start"])
    known = case.get("known") or []
    v1 = None
    if known:
        try:
            v1 = float(known[0]["viewers"])
        except (KeyError, TypeError, ValueError):
            v1 = None
    hist = _hourly_map(case.get("history_hourly") or [])
    if not hist:
        hist = _hourly_map(case.get("history") or [])
    if (not v1 or v1 <= 0) and hist:
        v1 = hist.get(max(hist))
    if not v1 or v1 <= 0:
        return None
    ln_v1 = math.log(v1)
    pre = [v for t, v in hist.items()
           if day - timedelta(days=7) <= t < day and v > 0]
    ln_prem = math.log(sum(pre) / len(pre)) if pre else ln_v1
    point12 = hist.get(day - timedelta(hours=12))
    g12 = None
    if point12 and point12 > 0:
        g12 = ln_v1 - math.log(point12)
    return {"ln_v1": ln_v1, "ln_prem": ln_prem, "g12": g12, "day": day}


def _build_library(case):
    lib = []
    for comp in case.get("comparables") or []:
        f = _launch_features(comp)
        if f is None:
            continue
        if f["tier"] == "target-history":
            continue
        if "grand theft auto vi" in f["name"].lower():
            continue
        lib.append(f)
    return lib


# ---------------------------------------------------------------------------
# Views: one OLS view (shared) and one-factor shape views per member
# ---------------------------------------------------------------------------

def _ols_view(lib, tgt):
    """Per-hour OLS view. Returns {hour: (median_log, residuals)}."""
    ln_v1s = [l["ln_v1"] for l in lib]
    lo, hi = min(ln_v1s), max(ln_v1s)
    clamp_v1 = min(max(tgt["ln_v1"], lo), hi + 0.5)
    g12s = [l["g12"] for l in lib if l["g12"] is not None]
    g12_med = _med(g12s) if g12s else 0.0
    tgt_g12 = tgt["g12"] if tgt["g12"] is not None else g12_med
    out = {}
    for h in HOURS:
        pool = []
        for l in lib:
            v = l["hourly"].get(l["day"] + timedelta(hours=h - 1))
            if v and v > 0:
                pool.append(l)
        if not pool:
            continue
        windowed = [l for l in pool if abs(l["ln_v1"] - clamp_v1) <= 2.5]
        use = windowed if len(windowed) >= 10 else pool
        if len(use) < 8:
            out[h] = (tgt["ln_v1"], None)
            continue
        rows = []
        ys = []
        for l in use:
            v = l["hourly"].get(l["day"] + timedelta(hours=h - 1))
            g = l["g12"] if l["g12"] is not None else g12_med
            rows.append([1.0, l["ln_v1"], l["ln_prem"], g])
            ys.append(math.log(v))
        beta = _lstsq(rows, ys)
        if beta is None:
            out[h] = (tgt["ln_v1"], None)
            continue
        pred = [sum(b * x for b, x in zip(beta, r)) for r in rows]
        resid = [y - p for y, p in zip(ys, pred)]
        centre = _med(resid)
        med = (beta[0] + beta[1] * clamp_v1 + beta[2] * tgt["ln_prem"]
               + beta[3] * tgt_g12 + centre)
        out[h] = (med, resid)
    return out, clamp_v1, g12_med


def _shape_view(pool, tgt, feats, sd_fn, windowed, clamp_v1):
    """One one-factor shape view. Returns {hour: (median_log, spread)}.

    feats: which target features drive the (level, amplitude) regressions.
    windowed: restrict the pool to +-2.5 log units around the target.
    """
    if windowed:
        pool = [l for l in pool if abs(l["ln_v1"] - clamp_v1) <= 2.5]
    if not pool:
        return {}
    template = {}
    for h in HOURS:
        vals = [l["ratios"][h] for l in pool if h in l["ratios"]]
        if vals:
            template[h] = _med(vals)
    if not template:
        return {}
    fits = {}
    for l in pool:
        xs = []
        ys = []
        for h, m in template.items():
            if h in l["ratios"]:
                xs.append(m)
                ys.append(l["ratios"][h])
        if len(xs) < 20:
            continue
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        den = sum((x - mx) ** 2 for x in xs)
        if den <= 1e-12:
            continue
        s = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
        c = my - s * mx
        fits[id(l)] = (c, s, l)
    feat_of = {
        "lnV1": lambda l: l["ln_v1"],
        "lnPREM": lambda l: l["ln_prem"],
        "g12": lambda l: (l["g12"] if l["g12"] is not None else 0.0),
    }
    tgt_g12 = tgt["g12"] if tgt["g12"] is not None else 0.0
    tgt_feats = {"lnV1": clamp_v1, "lnPREM": tgt["ln_prem"], "g12": tgt_g12}
    cs = None
    if len(fits) >= 5:
        rows_c = []
        rows_s = []
        ys_c = []
        ys_s = []
        for (c, s, l) in fits.values():
            g = l["g12"] if l["g12"] is not None else 0.0
            fvals = {"lnV1": l["ln_v1"], "lnPREM": l["ln_prem"], "g12": g}
            row = [1.0] + [fvals[f] for f in feats]
            rows_c.append(row)
            rows_s.append(row)
            ys_c.append(c)
            ys_s.append(s)
        beta_c = _lstsq(rows_c, ys_c)
        beta_s = _lstsq(rows_s, ys_s)
        if beta_c is not None and beta_s is not None:
            trow = [1.0] + [tgt_feats[f] for f in feats]
            c_star = sum(b * x for b, x in zip(beta_c, trow))
            s_star = sum(b * x for b, x in zip(beta_s, trow))
            cs = (c_star, min(max(s_star, 0.0), 2.0))
    if cs is None:
        all_c = [c for (c, s, l) in fits.values()]
        all_s = [s for (c, s, l) in fits.values()]
        if all_c:
            cs = (_med(all_c), min(max(_med(all_s), 0.0), 2.0))
        else:
            cs = (0.0, 1.0)
    c_star, s_star = cs
    out = {}
    for h, m in template.items():
        med = tgt["ln_v1"] + c_star + s_star * m
        resid = []
        for (c, s, l) in fits.values():
            if h in l["ratios"]:
                resid.append(l["ratios"][h] - (c + s * m))
        spread = sd_fn(resid) if len(resid) >= 2 else 0.4
        out[h] = (med, spread, m)
    return out


def _ea_keep(launch, kind):
    if kind == "tier":
        return launch["tier"] != "ea"
    if kind == "naive":
        return "early access" not in launch["unlock_low"]
    if kind == "hyphen":
        return "early-access" not in launch["unlock_low"]
    return True


def _member_values(lib, tgt, ols, clamp_v1, kind, shape_feats, sd_fn,
                   mix, spread_rule, windowed_shape):
    """One member's per-hour (log median, spread)."""
    if kind in ("tier", "naive"):
        shape_pool = [l for l in lib if _ea_keep(l, kind)]
    else:
        shape_pool = [l for l in lib if _ea_keep(l, "hyphen")]
    con_pool = [l for l in shape_pool if l["console"]]
    con = _shape_view(con_pool, tgt, shape_feats, sd_fn, windowed_shape,
                      clamp_v1)
    allv = _shape_view(shape_pool, tgt, shape_feats, sd_fn, windowed_shape,
                       clamp_v1)
    out = {}
    w_ols, w_con, w_all = mix
    for h in HOURS:
        if h not in ols:
            continue
        o_med, resid = ols[h]
        if resid is None:
            out[h] = (o_med, 0.4)
            continue
        s_ols = sd_fn(resid) * 0.8
        if h not in con or h not in allv:
            out[h] = (o_med, s_ols)
            continue
        c_med, c_sig, _ = con[h]
        a_med, a_sig, _ = allv[h]
        med = w_ols * o_med + w_con * c_med + w_all * a_med
        if spread_rule == "sqrt":
            var = (w_ols * s_ols ** 2 + w_con * c_sig ** 2
                   + w_all * a_sig ** 2)
            spread = max(math.sqrt(var), math.sqrt(0.02))
        else:
            spread = 0.81 * (w_ols * s_ols + w_con * c_sig + w_all * a_sig)
        out[h] = (med, spread)
    return out


def _member_values_y1(lib, tgt, ols, clamp_v1, shape_views):
    """Y1: same pools/spreads as Z0, but its shape (level, amplitude)
    regressions run on [lnV1, lnPREM] -- the v1 original features --
    instead of [lnV1, g12]."""
    con, allv = shape_views
    out = {}
    for h in HOURS:
        if h not in ols:
            continue
        o_med, resid = ols[h]
        if resid is None:
            out[h] = (o_med, 0.4)
            continue
        s_ols = _pstdev(resid) * 0.8
        if h not in con or h not in allv:
            out[h] = (o_med, s_ols)
            continue
        med = 0.5 * o_med + 0.2 * con[h][0] + 0.3 * allv[h][0]
        var = (0.5 * s_ols ** 2 + 0.2 * con[h][1] ** 2
               + 0.3 * allv[h][1] ** 2)
        spread = max(math.sqrt(var), math.sqrt(0.02))
        out[h] = (med, spread)
    return out


def members_for_case(case):
    """All four members for one case. Returns (members, ln_v1)."""
    tgt = _target_features(case)
    lib = _build_library(case)
    if tgt is None or not lib:
        return {}, None
    ols, clamp_v1, _ = _ols_view(lib, tgt)
    y5 = _member_values(lib, tgt, ols, clamp_v1, "tier",
                        ("lnV1", "g12"), _pstdev, (0.5, 0.2, 0.3),
                        "sqrt", True)
    z0 = _member_values(lib, tgt, ols, clamp_v1, "naive",
                        ("lnV1", "g12"), _pstdev, (0.5, 0.2, 0.3),
                        "sqrt", True)
    naive_pool = [l for l in lib if _ea_keep(l, "naive")]
    con_pool = [l for l in naive_pool if l["console"]]
    con_z = _shape_view(con_pool, tgt, ("lnV1", "lnPREM"), _pstdev, True,
                        clamp_v1)
    all_z = _shape_view(naive_pool, tgt, ("lnV1", "lnPREM"), _pstdev, True,
                        clamp_v1)
    y1 = _member_values_y1(lib, tgt, ols, clamp_v1, (con_z, all_z))
    clone = _member_values(lib, tgt, ols, clamp_v1, "hyphen",
                           ("lnV1", "lnPREM", "g12"), _stdev,
                           (0.6, 0.2, 0.2), "mean", False)
    members = {"y5": y5, "z0": z0, "y1": y1, "clone": clone}
    return members, tgt["ln_v1"]


# ---------------------------------------------------------------------------
# Combination: the packaged cell
# ---------------------------------------------------------------------------

def combine(members, ln_v1, params=None):
    """Per-hour (log median, lower spread, upper spread)."""
    prm = dict(P)
    if params:
        prm.update(params)
    out = {}
    for h in HOURS:
        vals = {k: members[k][h] for k in ("y5", "z0", "y1", "clone")
                if h in members.get(k, {})}
        if len(vals) < 4:
            continue
        m_y5, m_z0, m_y1, m_cl = (vals["y5"][0], vals["z0"][0],
                                 vals["y1"][0], vals["clone"][0])
        s_y5, s_z0, s_cl = vals["y5"][1], vals["z0"][1], vals["clone"][1]
        band = _band(h)
        if band == "early":
            if prm["early_src"] == "med3":
                mu0 = _med([m_y5, m_z0, m_y1])
            elif prm["early_src"] == "y1z0":
                mu0 = 0.5 * m_y1 + 0.5 * m_z0
            else:
                mu0 = m_y1
            ramp = _ramp_at(h, 0.0, prm["early_ramp"], 2, 6)
        elif band == "mid":
            if prm["mid_src"] == "med4":
                mu0 = _med([m_y5, m_z0, m_y1, m_cl])
            elif prm["mid_src"] == "y5clonez0":
                mu0 = 0.4 * m_y5 + 0.4 * m_cl + 0.2 * m_z0
            else:
                mu0 = 0.5 * m_y5 + 0.5 * m_cl
            ramp = _ramp_at(h, 0.0, prm["mid_ramp"], 7, 14)
        else:
            if prm["late_src"] == "med4":
                mu0 = _med([m_y5, m_z0, m_y1, m_cl])
            elif prm["late_src"] == "mean3":
                mu0 = (m_y5 + m_z0 + m_y1) / 3.0
            else:
                mu0 = _med([m_y5, m_z0, m_y1])
            ramp = _ramp_at(h, 0.0, prm["late_ramp"], 15, 24)
        size_x = max(0.0, ln_v1 - prm["size_pivot"])
        size_shift = prm["size_slope"] * size_x
        sd4 = _pstdev([m_y5, m_z0, m_y1, m_cl])
        sb = prm["shift_base"]
        if band == "late" and prm.get("shift_late") is not None:
            sb = prm["shift_late"]
        shift = (sb + prm["shift_slope"] * size_x) * sd4
        mu = mu0 + ramp + size_shift + shift
        w_y5, w_z0, w_cl = prm["spread_w"]
        if band == "late":
            base = _ramp_at(h, prm["spread_base"], prm["spread_late_end"],
                            15, 24)
        else:
            base = prm["spread_base"]
        spread = (base * (w_y5 * s_y5 + w_z0 * s_z0 + w_cl * s_cl)
                  * (1.0 + prm["disag_k"] * sd4)
                  * (1.0 + prm["spread_size_slope"] * size_x))
        lo_scale, hi_scale = prm["asym"][band]
        s_lo = min(spread * lo_scale, prm["spread_cap"])
        s_hi = min(spread * hi_scale, prm["spread_cap"])
        out[h] = (mu, s_lo, s_hi)
    return out


def _quantiles(mu, s_lo, s_hi):
    qs = [
        math.exp(mu - Z05 * s_lo),
        math.exp(mu - Z25 * s_lo),
        math.exp(mu),
        math.exp(mu + Z25 * s_hi),
        math.exp(mu + Z05 * s_hi),
    ]
    fixed = []
    prev = 0.0
    for q in qs:
        q = min(max(q, 1.0), 1e9)
        q = max(q, prev)
        fixed.append(q)
        prev = q
    return fixed


def forecast_case(case):
    members, ln_v1 = members_for_case(case)
    if not members or ln_v1 is None:
        raise ValueError("cannot build forecast for case")
    combined = combine(members, ln_v1)
    rows = []
    for h in HOURS:
        if h not in combined:
            raise ValueError(f"missing hour {h}")
        mu, s_lo, s_hi = combined[h]
        rows.append({"hour": h, "q": _quantiles(mu, s_lo, s_hi)})
    return {"hours": rows}


# ---------------------------------------------------------------------------
# Program entry point
# ---------------------------------------------------------------------------

def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR") or "/input"
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR") or "/output"
    os.makedirs(out_dir, exist_ok=True)
    names = sorted(n for n in os.listdir(in_dir) if n.endswith(".in"))
    if not names:
        return 0
    for name in names:
        with open(os.path.join(in_dir, name), "r", encoding="utf-8") as f:
            case = json.load(f)
        result = forecast_case(case)
        out_name = name[:-3] + ".out"
        with open(os.path.join(out_dir, out_name), "w",
                  encoding="utf-8") as f:
            json.dump(result, f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
