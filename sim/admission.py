#!/usr/bin/env python3
"""AgNP-channel admission ledger for a protocol-constant tip-ODE forcing.

Thesis #8 named the papaya AgNP assay as observation channels and barred
those scores from Theta. Thesis #19 admits some screen scores only as
forcing provenance. This script applies an AgNP-specific predicate to the
transcribed Thesis #8 channel records.

Kinetic Theta is the six decimal strings of Thesis #19. Reusing them is the
control: the digest must stay that digest. The carrier field is the same
demonstration field, reimplemented here so the separation can be checked.
Fisher matrices are not recomputed. Thesis 0 percentages are not in the
candidate file. No dose is selected.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
CANDIDATE_PATH = ROOT / "candidates.json"
# SHA-256 of sim/candidates.json as stored. A byte edit raises.
CANDIDATE_PIN = "6d7bb8b9d99302dd8d90361bc3d95f792b19ac0ddf99445693557dc97445ca5e"

FORCING_SYMBOL = "u_tip"
PROTOCOL_AMPLITUDE = "1"
HORIZON = "40"
THETA_NAMES = ("d0", "d_u", "g", "k_glyc", "s", "u_in")
# Published kinetic digest of Thesis #19. This ledger must reproduce it.
T19_THETA_SHA256 = "1ace2bb9eed19f85046f25ab690076d0d5920cf19219ddb362c775be9f41d0e8"

THETA = {
    "d0": "0.15",
    "d_u": "0.25",
    "g": "0.20",
    "k_glyc": "0.60",
    "s": "0.40",
    "u_in": "0.50",
}

G0, R0, A0 = 1.0, 0.20, 0.80
T_END = 40.0
N_SAMPLES = 401

AFFINE_ALPHA0 = 0.20
AFFINE_ALPHA1 = 0.50
PRIOR_SHIFT = 0.02
PRIOR_SD = 0.05
OBS_SD = 0.05

NUMERAL_KEYS = ("numeral", "rank", "leading_eigenvalue", "condition", "kappa")


def canon_bytes(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def digest(obj) -> str:
    return hashlib.sha256(canon_bytes(obj)).hexdigest()


def load_candidates() -> dict:
    raw = CANDIDATE_PATH.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != CANDIDATE_PIN:
        raise SystemExit(f"candidate pin mismatch: {got}")
    return json.loads(raw.decode("utf-8"))


def theta_float() -> dict:
    return {k: float(THETA[k]) for k in THETA_NAMES}


def record_by_id(doc: dict) -> dict:
    return {row["id"]: row for row in doc["records"]}


def numeral_strings(doc: dict) -> list[str]:
    found = []
    for row in doc["records"]:
        for key in NUMERAL_KEYS:
            value = row.get(key)
            if isinstance(value, str) and value:
                found.append(value)
    return found


def null_forcing() -> dict:
    return {
        "amplitude_rule": "null",
        "pieces": [{"t0": "0", "t1": HORIZON, "value": "0"}],
        "provenance": None,
        "symbol": FORCING_SYMBOL,
    }


def legal_forcing(rec: dict, results_sha: str) -> dict:
    return {
        "amplitude_rule": "protocol_constant",
        "pieces": [{"t0": "0", "t1": HORIZON, "value": PROTOCOL_AMPLITUDE}],
        "provenance": {
            "channel": rec["channel"],
            "class_label": "B",
            "named_agnp_channel": True,
            "results_sha256": results_sha,
            "score_copied": False,
            "score_part": "summary",
        },
        "symbol": FORCING_SYMBOL,
    }


def failed_checks(doc: dict, call: dict, slot_full: bool) -> list[str]:
    """Return every failed predicate, in the order Table 3-2 records them."""
    reasons: list[str] = []
    kind = call["kind"]
    destination = call["destination"]
    rule = call["amplitude_rule"]

    if kind != "channel_score":
        reasons.append("not_a_channel_score")
    if rule == "soft_prior" or kind == "soft_prior":
        reasons.append("soft_prior")
    if rule == "soft_weight" or kind == "soft_weight":
        reasons.append("soft_weight")
    if destination in THETA_NAMES:
        reasons.append("theta_destination")
    elif destination != FORCING_SYMBOL:
        reasons.append("symbol")

    if kind != "channel_score":
        if rule != "protocol_constant":
            reasons.append("amplitude_rule")
        if not reasons and slot_full:
            reasons.append("slot_occupied")
        return reasons

    rec = record_by_id(doc)[call["record_id"]]
    source = call.get("results_sha256", doc["results_sha256"])
    may = call.get("may_enter_theta", doc["may_enter_theta"])
    frozen = call.get("theta_frozen", doc["theta_frozen"])
    klass = call.get("class_label", rec.get("class_label", doc["class_label"]))
    free = call["free_theta"] if "free_theta" in call else rec["free_theta"]
    fitted = call.get("fitted_to_thesis0", doc["fitted_to_thesis0"])
    kappa_out = call.get("kappa_outside_theta", doc["kappa_outside_theta"])
    if "channel" in call and call["channel"] != rec["channel"]:
        reasons.append("channel_mismatch")
    if "score_part" in call and call["score_part"] != rec["score_part"]:
        reasons.append("part_mismatch")
    channel = rec["channel"]
    part = rec["score_part"]

    if source != doc["results_sha256"]:
        reasons.append("source_pin")
    if may is not False:
        reasons.append("may_enter_theta")
    if frozen is not True or klass != "B" or free is not False:
        reasons.append("class_A")
    if fitted is not False:
        reasons.append("wetlab_fit")
    if kappa_out is not True:
        reasons.append("kappa_in_theta")
    if channel != doc["named_agnp_channel"]:
        reasons.append("not_named_agnp_channel")
    if part != "summary":
        reasons.append("numeral_not_provenance")
    if rule != "protocol_constant":
        reasons.append("amplitude_rule")
    if not reasons and slot_full:
        reasons.append("slot_occupied")
    return reasons


def apply_call(doc: dict, ledger: dict, call: dict) -> dict:
    slot_full = ledger["forcing"]["amplitude_rule"] != "null"
    reasons = failed_checks(doc, call, slot_full=slot_full)
    record = {
        "id": call["id"],
        "kind": call["kind"],
        "record_id": call.get("record_id"),
        "destination": call["destination"],
        "amplitude_rule": call["amplitude_rule"],
        "reasons": reasons,
        "status": "refused" if reasons else "admitted",
        "theta_sha256": digest(ledger["theta"]),
        "note": call.get("note", ""),
    }
    if record["status"] == "admitted":
        rec = record_by_id(doc)[call["record_id"]]
        ledger["forcing"] = legal_forcing(rec, doc["results_sha256"])
        record["admitted_record"] = rec["id"]
    record["theta_sha256_after"] = digest(ledger["theta"])
    record["forcing_sha256_after"] = digest(ledger["forcing"])
    record["theta_unchanged"] = record["theta_sha256"] == record["theta_sha256_after"]
    ledger["calls"].append(record)
    return record


def equilibrium(theta: dict, u: float) -> dict:
    g_star = theta["u_in"] / theta["k_glyc"]
    r_star = theta["g"] * theta["u_in"] / (theta["d0"] + theta["d_u"] * u)
    a_star = theta["u_in"] / (theta["s"] * (1.0 + r_star))
    return {"G": g_star, "R": r_star, "A": a_star}


def closed_form(t: np.ndarray, theta: dict, u: float) -> tuple[np.ndarray, np.ndarray]:
    rate = theta["k_glyc"]
    g_inf = theta["u_in"] / rate
    g = g_inf + (G0 - g_inf) * np.exp(-rate * t)
    source = theta["g"] * rate
    clearance = theta["d0"] + theta["d_u"] * u
    if abs(clearance - rate) < 1e-12:
        raise SystemExit("closed form assumes the stress clearance differs from k_glyc")
    coef = source * (G0 - g_inf) / (clearance - rate)
    particular = source * g_inf / clearance
    homogeneous = R0 - coef - particular
    r = coef * np.exp(-rate * t) + particular + homogeneous * np.exp(-clearance * t)
    return g, r


def integrate(theta: dict, u: float) -> dict:
    t = np.linspace(0.0, T_END, N_SAMPLES)

    def field(_t, y):
        g, r, a_state = y
        d_g = theta["u_in"] - theta["k_glyc"] * g
        d_r = theta["g"] * theta["k_glyc"] * g - (theta["d0"] + theta["d_u"] * u) * r
        d_a = theta["k_glyc"] * g / (1.0 + r) - theta["s"] * a_state
        return (d_g, d_r, d_a)

    sol = solve_ivp(
        field,
        (0.0, T_END),
        (G0, R0, A0),
        t_eval=t,
        method="LSODA",
        rtol=1e-8,
        atol=1e-11,
    )
    if not sol.success:
        raise SystemExit(f"integration failed: {sol.message}")
    g_hat, r_hat = closed_form(sol.t, theta, u)
    g_err = float(np.max(np.abs(sol.y[0] - g_hat)))
    r_err = float(np.max(np.abs(sol.y[1] - r_hat)))
    if g_err > 1e-7 or r_err > 1e-7:
        raise SystemExit(f"closed form gap G {g_err} R {r_err}")
    if float(np.min(sol.y[1])) <= 0.0:
        raise SystemExit("stress coordinate left the positive half-line")
    packed = [
        [
            round(float(sol.t[i]), 10),
            round(float(sol.y[0, i]), 12),
            round(float(sol.y[1, i]), 12),
            round(float(sol.y[2, i]), 12),
        ]
        for i in range(sol.t.size)
    ]
    eq = equilibrium(theta, u)
    terminal = {
        "G": float(sol.y[0, -1]),
        "R": float(sol.y[1, -1]),
        "A": float(sol.y[2, -1]),
    }
    eq_gap = {k: abs(terminal[k] - eq[k]) for k in ("G", "R", "A")}
    return {
        "t": sol.t,
        "y": sol.y,
        "state_sha256": digest(packed),
        "terminal": terminal,
        "equilibrium": eq,
        "equilibrium_abs_gap": eq_gap,
        "closed_form_max_abs": {"G": g_err, "R": r_err},
        "u": u,
    }


def soft_prior_contrast(log_eigenvalue: float, d0: float) -> dict:
    """Normal-normal update in a side buffer. The ledger is not assigned."""
    mu = d0 + PRIOR_SHIFT * log_eigenvalue
    y_obs = d0
    post = (mu / PRIOR_SD**2 + y_obs / OBS_SD**2) / (1.0 / PRIOR_SD**2 + 1.0 / OBS_SD**2)
    leaked = dict(THETA)
    leaked["d0"] = f"{post:.16f}"
    return {
        "prior_mean": mu,
        "observation": y_obs,
        "posterior_mean_d0": post,
        "theta_sha256_if_written": digest(leaked),
        "written_to_ledger_theta": False,
    }


def eligibility(doc: dict) -> list[dict]:
    out = []
    for rec in doc["records"]:
        call = {
            "id": f"predicate-{rec['id']}",
            "kind": rec["kind"],
            "record_id": rec["id"],
            "destination": FORCING_SYMBOL,
            "amplitude_rule": "protocol_constant",
        }
        reasons = failed_checks(doc, call, slot_full=False)
        out.append(
            {
                "id": rec["id"],
                "kind": rec["kind"],
                "channel": rec.get("channel"),
                "score_part": rec.get("score_part"),
                "rank": rec.get("rank"),
                "eligible": not reasons,
                "reasons": reasons,
            }
        )
    return out


def build_calls(doc: dict) -> list[dict]:
    pin = doc["results_sha256"]
    broken = pin[:-1] + ("0" if pin[-1] != "0" else "1")
    symbol = FORCING_SYMBOL
    return [
        {
            "id": "R01",
            "kind": "theta_split",
            "record_id": "S_theta",
            "destination": "d0",
            "amplitude_rule": "copy_score",
            "note": "Class A inhibition 0.25 asked to replace d0",
        },
        {
            "id": "R02",
            "kind": "channel_score",
            "record_id": "S_lam",
            "destination": symbol,
            "amplitude_rule": "copy_eigenvalue",
            "note": "C_amyl leading eigenvalue asked to be the amplitude",
        },
        {
            "id": "R03",
            "kind": "channel_score",
            "record_id": "S_rank",
            "destination": "d0",
            "amplitude_rule": "soft_prior",
            "note": "Gaussian prior on d0 shifted by the C_amyl rank",
        },
        {
            "id": "R04",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "copy_eigenvalue",
            "note": "The eligible summary, with its eigenvalue copied into u_tip",
        },
        {
            "id": "R05",
            "kind": "channel_score",
            "record_id": "S_kc",
            "destination": symbol,
            "amplitude_rule": "copy_kappa",
            "note": "Declared K_c of C_amyl asked to be the amplitude",
        },
        {
            "id": "R06",
            "kind": "channel_score",
            "record_id": "S_const",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "note": "Constant channel asked to authorize the protocol schedule",
        },
        {
            "id": "R07",
            "kind": "channel_score",
            "record_id": "S_gate",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "note": "Gate channel, not the named AgNP channel",
        },
        {
            "id": "R08",
            "kind": "wetlab_percent",
            "record_id": "S_wet",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "note": "A wet-lab percent that this file does not contain",
        },
        {
            "id": "R09",
            "kind": "foreign_screen",
            "record_id": "S_l10",
            "destination": "u_phyto",
            "amplitude_rule": "protocol_constant",
            "note": "Thesis #19 row L10 aimed at Thesis #19's symbol",
        },
        {
            "id": "R10",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": "d0",
            "amplitude_rule": "protocol_constant",
            "note": "Legal-looking summary aimed at a kinetic name",
        },
        {
            "id": "R11",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": "u_phyto",
            "amplitude_rule": "protocol_constant",
            "note": "Legal summary aimed at a symbol this ledger does not declare",
        },
        {
            "id": "R12",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": "N_ROS",
            "amplitude_rule": "protocol_constant",
            "note": "Legal summary aimed at the input name of Thesis #7",
        },
        {
            "id": "R13",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "may_enter_theta": True,
            "note": "may_enter_theta flipped",
        },
        {
            "id": "R14",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "class_label": "A",
            "note": "Class A label written over the frozen class B record",
        },
        {
            "id": "R15",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "results_sha256": broken,
            "note": "One hex character of the Thesis #8 results pin flipped",
        },
        {
            "id": "R16",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "fitted_to_thesis0": True,
            "note": "Flag claiming the channel was fitted to Thesis 0",
        },
        {
            "id": "R17",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "kappa_outside_theta": False,
            "note": "Channel constant moved into Theta by a flag",
        },
        {
            "id": "R18",
            "kind": "channel_score",
            "record_id": "S_lam",
            "destination": "s",
            "amplitude_rule": "soft_weight",
            "note": "Eigenvalue used as a weight on the ATP clearance",
        },
        {
            "id": "A01",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "note": "Single legal admission of the C_amyl summary as provenance",
        },
        {
            "id": "R19",
            "kind": "channel_score",
            "record_id": "S_amyl",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "note": "Second write of the same legal record",
        },
        {
            "id": "R20",
            "kind": "channel_score",
            "record_id": "S_GA",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "note": "Joint channel after the slot is full",
        },
        {
            "id": "R21",
            "kind": "channel_score",
            "record_id": "S_const",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "channel": "C_amyl",
            "note": "Constant channel relabelled C_amyl inside the call",
        },
        {
            "id": "R22",
            "kind": "channel_score",
            "record_id": "S_lam",
            "destination": symbol,
            "amplitude_rule": "protocol_constant",
            "score_part": "summary",
            "note": "Eigenvalue record relabelled as a summary inside the call",
        },
    ]


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def draw_figures(elig: list[dict], calls: list[dict], curves: dict, eq_curve: dict) -> None:
    fig_dir = ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "figure.facecolor": "white",
        }
    )

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    labels = []
    for row in elig:
        channel = row["channel"] or row["kind"]
        labels.append(f"{row['id']}  {channel}  {row['score_part']}")
    colors = ["#1f4e79" if row["eligible"] else "#b9b9b9" for row in elig]
    y = np.arange(len(elig))
    ax.barh(y, [1 if row["eligible"] else 0 for row in elig], color=colors, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("AgNP admission predicate (1 = eligible, 0 = refused)")
    ax.set_xlim(0, 1.15)
    ax.invert_yaxis()
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / "eligibility.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.3), sharex=True)
    names = [("G", 0, "glucose-like"), ("R", 1, "stress-like"), ("A", 2, "ATP-like")]
    series = [
        ("null", curves["null"], "#555555", "-"),
        ("admitted u_tip = 1", curves["legal"], "#1f4e79", "-"),
        ("K_c copy, not admitted", curves["contrast"], "#a33c3c", ":"),
    ]
    for ax, (key, idx, title) in zip(axes, names):
        for name, curve, color, ls in series:
            ax.plot(curve["t"], curve["y"][idx], color=color, ls=ls, lw=1.6, label=name)
        ax.set_title(title)
        ax.set_xlabel("toy time")
        ax.set_ylabel(key)
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "trajectories.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    xs = np.arange(len(calls))
    ax.scatter(xs, np.ones(len(calls)), s=28, color="#1f4e79", label="Θ digest unchanged", zorder=3)
    refused_x = [i for i, c in enumerate(calls) if c["status"] != "admitted"]
    admitted_x = [i for i, c in enumerate(calls) if c["status"] == "admitted"]
    ax.scatter(
        refused_x,
        np.zeros(len(refused_x)),
        s=22,
        color="#9a9a9a",
        label="forcing not written by this call",
        zorder=2,
    )
    if admitted_x:
        ax.scatter(
            admitted_x,
            np.zeros(len(admitted_x)),
            s=42,
            color="#a33c3c",
            label="forcing schedule written",
            zorder=4,
        )
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["forcing slot", "Θ digest"])
    ax.set_xticks(xs)
    ax.set_xticklabels([c["id"] for c in calls], rotation=70, fontsize=7)
    ax.set_ylim(-0.25, 1.35)
    ax.set_xlabel("ledger call")
    style_axes(ax)
    ax.legend(frameon=False, fontsize=7, loc="center right")
    fig.tight_layout()
    fig.savefig(fig_dir / "digest_stability.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(eq_curve["u"], eq_curve["R"], color="#1f4e79", lw=1.6, label="R equilibrium")
    ax.plot(eq_curve["u"], eq_curve["A"], color="#3d6b4f", lw=1.6, label="A equilibrium")
    ax.axvline(1.0, color="#1f4e79", ls="--", lw=0.9, label="protocol amplitude 1")
    ax.scatter(
        [eq_curve["copy_u"]],
        [eq_curve["copy_R"]],
        color="#a33c3c",
        s=28,
        zorder=3,
        label="K_c copy, not admitted",
    )
    ax.set_xlabel("constant u_tip")
    ax.set_ylabel("equilibrium")
    style_axes(ax)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "equilibria.png", dpi=140)
    plt.close(fig)


def snapshot_times(curve: dict) -> dict:
    want = (0.0, 5.0, 10.0, 20.0, 40.0)
    out = {}
    for w in want:
        i = int(np.argmin(np.abs(curve["t"] - w)))
        out[f"{w:.0f}"] = {
            "G": round(float(curve["y"][0, i]), 6),
            "R": round(float(curve["y"][1, i]), 6),
            "A": round(float(curve["y"][2, i]), 6),
        }
    return out


def pack_curve(curve: dict) -> dict:
    return {
        "u": curve["u"],
        "state_sha256": curve["state_sha256"],
        "terminal": {k: round(v, 6) for k, v in curve["terminal"].items()},
        "equilibrium": {k: round(v, 6) for k, v in curve["equilibrium"].items()},
        "equilibrium_abs_gap": {k: float(f"{v:.6e}") for k, v in curve["equilibrium_abs_gap"].items()},
        "closed_form_max_abs": {k: float(f"{v:.6e}") for k, v in curve["closed_form_max_abs"].items()},
        "samples": snapshot_times(curve),
    }


def main() -> None:
    doc = load_candidates()
    if doc["results_sha256"] != "a62bd39023418c48c3b3fdff13cde90e6a136ea7a813ab0ab935feff6d985697":
        raise SystemExit("Thesis #8 results pin drifted")
    if doc["named_agnp_channel"] != "C_amyl":
        raise SystemExit("named channel drifted")
    if doc["wetlab_percent_transcribed"] is not False:
        raise SystemExit("a wet-lab percent was transcribed")
    if list(THETA) != list(THETA_NAMES):
        raise SystemExit("theta key order drifted")
    if digest(THETA) != T19_THETA_SHA256:
        raise SystemExit(f"theta digest is not the Thesis #19 digest: {digest(THETA)}")

    elig = eligibility(doc)
    eligible_ids = [row["id"] for row in elig if row["eligible"]]
    if eligible_ids != ["S_amyl"]:
        raise SystemExit(f"eligible set {eligible_ids}")

    summaries = [row for row in elig if row["score_part"] == "summary"]
    for row in summaries:
        if row["id"] == "S_amyl":
            continue
        if row["reasons"] != ["not_named_agnp_channel"]:
            raise SystemExit(f"companion reason drifted for {row['id']}: {row['reasons']}")
    for row in elig:
        if row["channel"] == "C_amyl" and row["score_part"] != "summary":
            if row["reasons"] != ["numeral_not_provenance"]:
                raise SystemExit(f"numeral reason drifted for {row['id']}: {row['reasons']}")

    ledger = {"theta": dict(THETA), "forcing": null_forcing(), "calls": []}
    theta_baseline = digest(ledger["theta"])
    forcing_null = digest(ledger["forcing"])
    theta_bytes = canon_bytes(ledger["theta"])

    for call in build_calls(doc):
        apply_call(doc, ledger, call)

    if canon_bytes(ledger["theta"]) != theta_bytes:
        raise SystemExit("theta bytes changed")
    if digest(ledger["theta"]) != theta_baseline:
        raise SystemExit("theta digest changed")
    admitted = [c for c in ledger["calls"] if c["status"] == "admitted"]
    if len(admitted) != 1 or admitted[0]["id"] != "A01":
        raise SystemExit(f"admission count {len(admitted)}")
    if any(not c["theta_unchanged"] for c in ledger["calls"]):
        raise SystemExit("a call rewrote theta")
    forcing_legal = digest(ledger["forcing"])
    if forcing_legal == forcing_null:
        raise SystemExit("legal forcing digest did not move")
    if ledger["forcing"]["provenance"]["score_copied"] is not False:
        raise SystemExit("provenance says the score was copied")
    if ledger["forcing"]["provenance"]["channel"] != "C_amyl":
        raise SystemExit("provenance channel is not C_amyl")

    payload = canon_bytes(ledger["forcing"]).decode("utf-8")
    for numeral in numeral_strings(doc):
        if numeral in payload:
            raise SystemExit(f"numeral {numeral} leaked into the forcing payload")

    pre = [c for c in ledger["calls"] if c["id"].startswith("R") and int(c["id"][1:]) <= 18]
    if any(c["forcing_sha256_after"] != forcing_null for c in pre):
        raise SystemExit("a refusal wrote the forcing")
    post = [c for c in ledger["calls"] if c["id"] in {"R19", "R20", "R21", "R22"}]
    if any(c["forcing_sha256_after"] != forcing_legal for c in post):
        raise SystemExit("a later refusal replaced the admitted schedule")
    by_call = {c["id"]: c for c in ledger["calls"]}
    if by_call["R19"]["reasons"] != ["slot_occupied"]:
        raise SystemExit(f"R19 reasons {by_call['R19']['reasons']}")
    if "slot_occupied" in by_call["R20"]["reasons"]:
        raise SystemExit("R20 was labelled only as occupancy")
    if by_call["R21"]["reasons"] != ["channel_mismatch", "not_named_agnp_channel"]:
        raise SystemExit(f"R21 reasons {by_call['R21']['reasons']}")
    if by_call["R22"]["reasons"] != ["part_mismatch", "numeral_not_provenance"]:
        raise SystemExit(f"R22 reasons {by_call['R22']['reasons']}")

    by_id = record_by_id(doc)
    amyl = by_id["S_amyl"]
    log_lam = math.log10(float(amyl["leading_eigenvalue"]))
    copy_u = float(by_id["S_kc"]["numeral"])
    affine_u = AFFINE_ALPHA0 + AFFINE_ALPHA1 * log_lam

    theta = theta_float()
    curve_null = integrate(theta, 0.0)
    curve_legal = integrate(theta, float(PROTOCOL_AMPLITUDE))
    curve_copy = integrate(theta, copy_u)
    curve_affine = integrate(theta, affine_u)
    if curve_null["state_sha256"] == curve_legal["state_sha256"]:
        raise SystemExit("state digest did not move under the admitted forcing")
    if digest(THETA) != theta_baseline:
        raise SystemExit("integration rewrote theta")

    prior = soft_prior_contrast(log_lam, float(THETA["d0"]))
    if prior["theta_sha256_if_written"] == theta_baseline:
        raise SystemExit("soft prior contrast would not have changed the digest")
    if prior["written_to_ledger_theta"]:
        raise SystemExit("contrast was written")
    if digest(ledger["theta"]) != theta_baseline:
        raise SystemExit("contrast touched the ledger")

    eq0 = equilibrium(theta, 0.0)
    eq1 = equilibrium(theta, 1.0)
    if abs(eq0["G"] - eq1["G"]) > 1e-15:
        raise SystemExit("forcing moved the glucose equilibrium")
    if abs(eq0["R"] - (0.20 * 0.50) / 0.15) > 1e-12:
        raise SystemExit("null stress equilibrium")
    if abs(eq1["R"] - 0.25) > 1e-12 or abs(eq1["A"] - 1.0) > 1e-12:
        raise SystemExit("admitted equilibrium drifted")

    max_state_gap = float(np.max(np.abs(curve_legal["y"] - curve_null["y"])))
    max_g_gap = float(np.max(np.abs(curve_legal["y"][0] - curve_null["y"][0])))

    u_grid = np.linspace(0.0, 2.5, 101)
    eq_curve = {
        "u": u_grid,
        "R": np.array([equilibrium(theta, float(u))["R"] for u in u_grid]),
        "A": np.array([equilibrium(theta, float(u))["A"] for u in u_grid]),
        "copy_u": copy_u,
        "copy_R": equilibrium(theta, copy_u)["R"],
    }
    draw_figures(
        elig,
        ledger["calls"],
        {"null": curve_null, "legal": curve_legal, "contrast": curve_copy},
        eq_curve,
    )

    refused = [c for c in ledger["calls"] if c["status"] == "refused"]
    results = {
        "rng_draws": 0,
        "candidate_sha256": CANDIDATE_PIN,
        "t08_results_sha256_transcribed": doc["results_sha256"],
        "t19_theta_sha256_control": T19_THETA_SHA256,
        "named_agnp_channel": "C_amyl",
        "theta": THETA,
        "theta_sha256": theta_baseline,
        "theta_sha256_after_all_calls": digest(ledger["theta"]),
        "theta_unchanged": theta_baseline == digest(ledger["theta"]) == T19_THETA_SHA256,
        "forcing_symbol": FORCING_SYMBOL,
        "forcing_null_sha256": forcing_null,
        "forcing_admitted_sha256": forcing_legal,
        "admitted_forcing": ledger["forcing"],
        "n_calls": len(ledger["calls"]),
        "n_refused": len(refused),
        "n_admitted": len(admitted),
        "eligible_ids": eligible_ids,
        "eligibility": elig,
        "calls": [
            {
                "id": c["id"],
                "kind": c["kind"],
                "record_id": c["record_id"],
                "destination": c["destination"],
                "amplitude_rule": c["amplitude_rule"],
                "status": c["status"],
                "reasons": c["reasons"],
                "theta_unchanged": c["theta_unchanged"],
                "note": c["note"],
            }
            for c in ledger["calls"]
        ],
        "curves": {
            "null": pack_curve(curve_null),
            "admitted": pack_curve(curve_legal),
            "kappa_copy_contrast_not_a_ledger_write": pack_curve(curve_copy),
            "affine_contrast_not_a_ledger_write": pack_curve(curve_affine),
        },
        "max_abs_state_gap_admitted_vs_null": round(max_state_gap, 6),
        "max_abs_G_gap_admitted_vs_null": float(f"{max_g_gap:.6e}"),
        "contrasts_not_written": {
            "soft_prior": {
                "record_id": "S_amyl",
                "column": "log10(leading_eigenvalue)",
                "log10_leading_eigenvalue": log_lam,
                "prior_mean_d0": prior["prior_mean"],
                "posterior_mean_d0": prior["posterior_mean_d0"],
                "theta_sha256_if_written": prior["theta_sha256_if_written"],
                "written_to_ledger_theta": False,
                "ledger_theta_sha256_after_contrast": digest(ledger["theta"]),
            },
            "kappa_copy_amplitude": {
                "record_id": "S_kc",
                "numeral": by_id["S_kc"]["numeral"],
                "proposed_u": copy_u,
                "written_to_ledger_forcing": False,
            },
            "affine_amplitude": {
                "record_id": "S_lam",
                "map": "0.20 + 0.50 * log10(leading_eigenvalue)",
                "proposed_u": affine_u,
                "written_to_ledger_forcing": False,
            },
        },
        "checks": {
            "eligible_set_is_S_amyl_only": True,
            "companion_summaries_fail_only_not_named_agnp_channel": True,
            "c_amyl_numerals_fail_only_numeral_not_provenance": True,
            "single_admission_A01": True,
            "theta_digest_equals_thesis_19": True,
            "theta_bytes_unchanged": True,
            "refusals_before_admission_left_null_forcing": True,
            "refusals_after_admission_left_admitted_forcing": True,
            "numerals_absent_from_forcing_payload": True,
            "closed_form_G_R": True,
            "G_equilibrium_independent_of_u": True,
            "state_digest_moves_when_forcing_is_admitted": True,
            "soft_prior_digest_differs_and_was_not_written": True,
            "wetlab_percent_not_transcribed": True,
        },
    }
    out = ROOT / "results.json"
    out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"theta_sha256 {theta_baseline}")
    print(f"forcing_null {forcing_null}")
    print(f"forcing_admitted {forcing_legal}")
    print(f"refused {len(refused)} admitted {len(admitted)} calls {len(ledger['calls'])}")
    print(f"eligible {eligible_ids}")
    print(f"state null {curve_null['state_sha256']}")
    print(f"state legal {curve_legal['state_sha256']}")
    print(f"max gap {max_state_gap:.6f} G gap {max_g_gap:.6e}")
    print(f"log10_lam {log_lam:.12f}")
    print(f"prior mu {prior['prior_mean']:.12f} post {prior['posterior_mean_d0']:.12f}")
    print(f"prior sha {prior['theta_sha256_if_written']}")
    print(f"copy u {copy_u:.6f} affine u {affine_u:.12f}")
    print(f"R19 {by_call['R19']['reasons']}")
    print(f"R21 {by_call['R21']['reasons']}")
    print(f"R22 {by_call['R22']['reasons']}")
    for c in ledger["calls"]:
        print(f"{c['id']} {c['status']} {c['reasons']}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
