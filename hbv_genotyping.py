"""HBV genotyping and sliding-window recombination screening.

The whole-genome classifier uses an Aurora-controlled HBV reference database.
A second windowed screen deliberately ignores genotype I and putative genotype J
when assigning parental-like segments, so a known recombinant lineage such as
HBV genotype I can still show its A/C/G-like mosaic signal.

This is a screening classifier, not a replacement for formal phylogenetic,
bootscan, jpHMM, or RDP confirmation of unusual recombinants.
"""
from __future__ import annotations

import csv
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from io import StringIO

from .cancellation import run_cancellable


HBV_GENOTYPE_DB_ENV = "AURORA_HBV_GENOTYPE_DB"
HBV_RECOMB_WINDOW_ENV = "AURORA_HBV_RECOMB_WINDOW"
HBV_RECOMB_STEP_ENV = "AURORA_HBV_RECOMB_STEP"
HBV_RECOMB_MIN_IDENTITY_ENV = "AURORA_HBV_RECOMB_MIN_IDENTITY"
HBV_RECOMB_MIN_QCOV_ENV = "AURORA_HBV_RECOMB_MIN_QCOV"
HBV_RECOMB_MIN_MARGIN_ENV = "AURORA_HBV_RECOMB_MIN_IDENTITY_MARGIN"
HBV_RECOMB_ALWAYS_ENV = "AURORA_HBV_RECOMB_ALWAYS"
HBV_RECOMB_TRIGGER_BIT_MARGIN_ENV = "AURORA_HBV_RECOMB_TRIGGER_BIT_MARGIN"
HBV_RECOMB_TRIGGER_IDENTITY_MARGIN_ENV = "AURORA_HBV_RECOMB_TRIGGER_IDENTITY_MARGIN"
HBV_RECOMB_MIN_COMPONENT_FRACTION_ENV = "AURORA_HBV_RECOMB_MIN_COMPONENT_FRACTION"
HBV_RECOMB_MIN_COMPONENT_RUN_ENV = "AURORA_HBV_RECOMB_MIN_COMPONENT_RUN"

DEFAULT_WINDOW = 400
DEFAULT_STEP = 100
DEFAULT_MIN_IDENTITY = 82.0
DEFAULT_MIN_QCOV = 80.0
DEFAULT_MIN_IDENTITY_MARGIN = 1.0
DEFAULT_TRIGGER_BIT_MARGIN = 1.0
DEFAULT_TRIGGER_IDENTITY_MARGIN = 1.0
DEFAULT_MIN_COMPONENT_FRACTION = 0.20
DEFAULT_MIN_COMPONENT_RUN = 4

# I is a recognized complex recombinant lineage. J remains putative and is
# excluded from the parental-like scan as well. Whole-genome genotyping still
# allows both when they are present in the HBV genotype reference database.
PARENTAL_SCAN_GENOTYPES = tuple("ABCDEFGH")
GENOTYPE_I_PARENTAL_SET = {"A", "C", "G"}


def _package_dir():
    return os.path.dirname(os.path.abspath(__file__))


def _db_prefix_exists(prefix):
    if not prefix:
        return False
    return any(os.path.exists(prefix + suffix) for suffix in (
        ".nal", ".ndb", ".nin", ".00.nin", ".nsq", ".00.nsq"
    ))


def find_hbv_genotype_database():
    candidates = []
    explicit = os.environ.get(HBV_GENOTYPE_DB_ENV, "").strip()
    if explicit:
        candidates.append(explicit)
    candidates.append(os.path.join(_package_dir(), "blastdb", "hbv_genotype_db", "hbv_genotypes"))
    seen = set()
    for candidate in candidates:
        if candidate and candidate not in seen:
            seen.add(candidate)
            if _db_prefix_exists(candidate):
                return candidate
    return ""


def _read_single_fasta_sequence(path):
    parts = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    parts.append(line.strip())
    except Exception:
        return ""
    return re.sub(r"[^ACGTN]", "N", "".join(parts).upper())


def _parse_label(title):
    text = str(title or "")
    g = re.search(r"\bgenotype=([A-J])\b", text, re.I)
    if not g:
        return "", "", ""
    sg = re.search(r"\bsubgenotype=([A-J][0-9A-Za-z.-]*)\b", text, re.I)
    acc = re.search(r"\baccession=([A-Z]{1,4}_?\d+(?:\.\d+)?)\b", text, re.I)
    return (
        g.group(1).upper(),
        sg.group(1).upper() if sg else "",
        acc.group(1).upper() if acc else "",
    )


def _parse_output(text):
    hits = []
    if not str(text or "").strip():
        return hits
    for fields in csv.reader(StringIO(text), delimiter="\t"):
        if len(fields) < 8:
            continue
        try:
            genotype, subgenotype, source_accession = _parse_label(fields[7])
            hits.append({
                "qseqid": fields[0],
                "subject_id": fields[1],
                "pident": float(fields[2]),
                "qcovs": float(fields[3]),
                "length": int(float(fields[4])),
                "evalue": fields[5],
                "bitscore": float(fields[6]),
                "stitle": fields[7],
                "genotype": genotype,
                "subgenotype": subgenotype,
                "source_accession": source_accession,
            })
        except (ValueError, IndexError):
            continue
    return hits


def _run_blast(blastn, query, db_prefix, threads, timeout, task="blastn", max_targets=100):
    outfmt = "6 qseqid sseqid pident qcovs length evalue bitscore stitle"
    cmd = [
        blastn,
        "-task", task,
        "-query", query,
        "-db", db_prefix,
        "-max_target_seqs", str(max_targets),
        "-evalue", "1e-10",
        "-num_threads", str(max(1, int(threads))),
        "-outfmt", outfmt,
    ]
    if task == "blastn":
        cmd.extend(["-dust", "no"])
    return run_cancellable(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )


def _best_per_genotype(hits, allowed=None):
    best = {}
    allowed_set = set(allowed) if allowed else None
    for hit in hits:
        genotype = hit.get("genotype", "")
        if not genotype or (allowed_set is not None and genotype not in allowed_set):
            continue
        old = best.get(genotype)
        score = (
            float(hit.get("bitscore", 0) or 0),
            float(hit.get("qcovs", 0) or 0),
            float(hit.get("pident", 0) or 0),
        )
        if old is None:
            best[genotype] = hit
        else:
            old_score = (
                float(old.get("bitscore", 0) or 0),
                float(old.get("qcovs", 0) or 0),
                float(old.get("pident", 0) or 0),
            )
            if score > old_score:
                best[genotype] = hit
    return best


def _rank_genotypes(best_by_genotype):
    return sorted(
        best_by_genotype.values(),
        key=lambda h: (
            float(h.get("bitscore", 0) or 0),
            float(h.get("qcovs", 0) or 0),
            float(h.get("pident", 0) or 0),
        ),
        reverse=True,
    )


def _margin(top, second):
    if not top or not second:
        return 999.0, 999.0
    top_bits = float(top.get("bitscore", 0) or 0)
    second_bits = float(second.get("bitscore", 0) or 0)
    bit_pct = ((top_bits - second_bits) / top_bits * 100.0) if top_bits else 0.0
    ident_pp = float(top.get("pident", 0) or 0) - float(second.get("pident", 0) or 0)
    return bit_pct, ident_pp


def _sequence_quality(sequence):
    seq = re.sub(r"[^ACGTN]", "", str(sequence or "").upper())
    total = len(seq)
    n_count = seq.count("N")
    callable_bases = total - n_count
    callable_percent = (callable_bases / total * 100.0) if total else 0.0
    gaps = [len(m.group(0)) for m in re.finditer(r"N+", seq)]
    return {
        "length": total,
        "n_count": n_count,
        "callable_bases": callable_bases,
        "callable_percent": callable_percent,
        "largest_n_gap": max(gaps) if gaps else 0,
    }


def _whole_genome_call(hits, sequence_quality=None):
    quality = sequence_quality or {}
    callable_bases = int(quality.get("callable_bases", 0) or 0)
    callable_percent = float(quality.get("callable_percent", 0.0) or 0.0)
    largest_n_gap = int(quality.get("largest_n_gap", 0) or 0)

    ranked = _rank_genotypes(_best_per_genotype(hits))
    if not ranked:
        return {
            "genotype": "", "subgenotype": "", "display": "Indeterminate",
            "confidence": "Low", "support": "Ingen genotypmärkt referensträff erhölls.",
            "best_hit": None, "second_hit": None,
            "bit_margin": 0.0, "identity_margin": 0.0,
        }

    top = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None
    bit_margin, ident_margin = _margin(top, second)
    qcov = float(top.get("qcovs", 0) or 0)
    ident = float(top.get("pident", 0) or 0)

    # Incomplete genomes must not automatically become E/D, B/C, etc. A partial
    # genome can still support one genotype when the available sequence is long
    # enough, highly similar, and clearly separated from the runner-up.
    insufficient_sequence = callable_bases < 1000
    weak_match = ident < 85.0
    near_tie = bool(second) and bit_margin < 0.5 and ident_margin < 0.5

    if insufficient_sequence or weak_match or near_tie:
        pair = "/".join(x.get("genotype", "?") for x in ranked[:2])
        reason = []
        if insufficient_sequence:
            reason.append(f"endast {callable_bases} analyserbara baser")
        if weak_match:
            reason.append(f"låg identitet {ident:.3f}%")
        if near_tie:
            reason.append(f"liten separation till nästa genotyp ({bit_margin:.2f}% bit-score-marginal)")
        support = "Genotypningen är osäker" + (": " + "; ".join(reason) if reason else ".")
        return {
            "genotype": "", "subgenotype": "", "display": f"Indeterminate ({pair})",
            "confidence": "Low", "support": support,
            "best_hit": top, "second_hit": second,
            "bit_margin": bit_margin, "identity_margin": ident_margin,
        }

    genotype = top.get("genotype", "")
    subgenotype = top.get("subgenotype", "")
    display = genotype if (not subgenotype or str(subgenotype).upper() == str(genotype).upper()) else f"{genotype} ({subgenotype})"

    clear_separation = (bit_margin >= 1.0 or ident_margin >= 0.75)
    if qcov >= 90.0 and callable_percent >= 90.0 and ident >= 90.0 and clear_separation:
        confidence = "High"
    elif qcov >= 80.0 and callable_percent >= 80.0 and ident >= 88.0 and clear_separation:
        confidence = "Moderate"
    else:
        # A genotype may still be reportable from a substantial, concordant
        # partial genome. This confidence explicitly prevents whole-genome or
        # recombination claims from being inferred from missing sequence.
        confidence = "Limited"

    support = (
        f"HBV-referensjämförelse: identitet {ident:.3f}%, BLAST query coverage {qcov:.1f}%, "
        f"analyserbar genomandel {callable_percent:.1f}% ({callable_bases} baser), "
        f"största N-gap {largest_n_gap} nt, bit-score-marginal {bit_margin:.2f}% till nästa genotyp."
    )
    return {
        "genotype": genotype,
        "subgenotype": subgenotype,
        "display": display,
        "confidence": confidence,
        "support": support,
        "best_hit": top,
        "second_hit": second,
        "bit_margin": bit_margin,
        "identity_margin": ident_margin,
    }


def _env_int(name, default, minimum=1):
    try:
        return max(minimum, int(os.environ.get(name, default)))
    except (TypeError, ValueError):
        return default


def _env_float(name, default):
    try:
        return float(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


def _env_bool(name, default=False):
    value = str(os.environ.get(name, "")).strip().lower()
    if not value:
        return bool(default)
    return value in {"1", "true", "yes", "on"}


def _recombination_indicated(whole):
    """Decide whether the expensive/windowed recombination screen is warranted.

    Genotype I is always characterized because recombinant ancestry is part of
    the lineage. Conventional A-H calls are screened only when the whole-genome
    result is indeterminate/moderate or when the best-vs-next genotype margin is
    unusually small. Set AURORA_HBV_RECOMB_ALWAYS=1 during validation to force
    the window scan for every HBV sample.
    """
    if _env_bool(HBV_RECOMB_ALWAYS_ENV, False):
        return True, "Forced by AURORA_HBV_RECOMB_ALWAYS."

    genotype = str(whole.get("genotype", "") or "").upper()
    confidence = str(whole.get("confidence", "") or "").lower()
    if genotype == "I":
        return True, "Genotype I is a recognized recombinant lineage and is always characterized."
    if not genotype:
        return True, "Whole-genome genotype is indeterminate."
    if confidence != "high":
        return True, f"Whole-genome genotype confidence is {whole.get('confidence', 'not high')}."

    bit_trigger = _env_float(HBV_RECOMB_TRIGGER_BIT_MARGIN_ENV, DEFAULT_TRIGGER_BIT_MARGIN)
    ident_trigger = _env_float(HBV_RECOMB_TRIGGER_IDENTITY_MARGIN_ENV, DEFAULT_TRIGGER_IDENTITY_MARGIN)
    bit_margin = float(whole.get("bit_margin", 0.0) or 0.0)
    ident_margin = float(whole.get("identity_margin", 0.0) or 0.0)

    # Require both margins to be narrow before a high-confidence conventional
    # genotype is sent into the recombination screen. This avoids triggering on
    # isolated locally similar regions in otherwise clear A-H genomes.
    if bit_margin < bit_trigger and ident_margin < ident_trigger:
        return True, (
            f"Whole-genome best-vs-next genotype separation is narrow "
            f"(bit-score margin {bit_margin:.2f}%, identity margin {ident_margin:.3f} pp)."
        )

    return False, (
        f"High-confidence genotype {genotype} with clear whole-genome separation "
        f"(bit-score margin {bit_margin:.2f}%, identity margin {ident_margin:.3f} pp)."
    )


def _make_windows(sequence, window, step):
    windows = []
    n = len(sequence)
    if n < max(250, window // 2):
        return windows
    start = 0
    idx = 1
    while start < n:
        end = min(start + window, n)
        piece = sequence[start:end]
        if len(piece) < max(250, int(window * 0.70)):
            break
        n_fraction = piece.count("N") / max(len(piece), 1)
        windows.append({
            "id": f"w{idx:03d}_{start + 1}_{end}",
            "start": start + 1,
            "end": end,
            "sequence": piece,
            "n_fraction": n_fraction,
        })
        idx += 1
        start += step
    return windows


def _assign_window(window, hits, min_identity, min_qcov, min_identity_margin):
    usable = [
        h for h in hits
        if h.get("genotype") in PARENTAL_SCAN_GENOTYPES
        and float(h.get("qcovs", 0) or 0) >= min_qcov
        and float(h.get("pident", 0) or 0) >= min_identity
    ]
    ranked = _rank_genotypes(_best_per_genotype(usable, PARENTAL_SCAN_GENOTYPES))
    if not ranked:
        return {**window, "genotype": "Ambiguous", "identity": "", "qcov": "", "margin": "", "accession": ""}
    top = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None
    bit_margin, ident_margin = _margin(top, second)

    # A window is assigned only when the leading genotype is clearly separated
    # from the runner-up. This is intentionally conservative to avoid false
    # recombinant calls from short, conserved HBV regions.
    clear = second is None or ident_margin >= min_identity_margin or bit_margin >= 1.5
    genotype = top.get("genotype", "") if clear else "Ambiguous"
    return {
        **window,
        "genotype": genotype,
        "identity": round(float(top.get("pident", 0) or 0), 3),
        "qcov": round(float(top.get("qcovs", 0) or 0), 1),
        "margin": round(float(ident_margin), 3) if second else "",
        "accession": top.get("source_accession", "") or top.get("subject_id", ""),
    }


def _component_stats(assignments):
    assigned = [x for x in assignments if x.get("genotype") in PARENTAL_SCAN_GENOTYPES]
    counts = Counter(x["genotype"] for x in assigned)
    longest_runs = Counter()

    previous = None
    run_len = 0
    for item in assignments + [{"genotype": None}]:
        genotype = item.get("genotype")
        if genotype in PARENTAL_SCAN_GENOTYPES and genotype == previous:
            run_len += 1
        else:
            if previous in PARENTAL_SCAN_GENOTYPES:
                longest_runs[previous] = max(longest_runs[previous], run_len)
            previous = genotype if genotype in PARENTAL_SCAN_GENOTYPES else None
            run_len = 1 if previous else 0
    return assigned, counts, longest_runs


def _strong_components(assignments, whole_genotype=""):
    assigned, counts, longest_runs = _component_stats(assignments)
    if not assigned:
        return [], counts, longest_runs

    min_fraction = _env_float(
        HBV_RECOMB_MIN_COMPONENT_FRACTION_ENV, DEFAULT_MIN_COMPONENT_FRACTION
    )
    min_run = _env_int(
        HBV_RECOMB_MIN_COMPONENT_RUN_ENV, DEFAULT_MIN_COMPONENT_RUN, minimum=2
    )
    min_count = max(min_run, int(math.ceil(len(assigned) * min_fraction)))

    strong = [
        g for g, count in counts.most_common()
        if count >= min_count and longest_runs.get(g, 0) >= min_run
    ]
    return strong, counts, longest_runs


def _classify_recombination(assignments, whole_genotype):
    if not assignments:
        return "Indeterminate", "Windowed recombination screen could not be performed.", ""

    strong, counts, longest_runs = _strong_components(assignments, whole_genotype)
    assigned_count = sum(counts.values())
    ambiguous_count = sum(1 for x in assignments if x.get("genotype") == "Ambiguous")
    components_text = "/".join(strong) if strong else ""

    if whole_genotype == "I":
        # Genotype I is a recognized recombinant lineage, so a nearest-reference
        # whole-genome hit is treated as a *candidate*, not sufficient proof on
        # its own.  Require the independent parental-window screen to recover at
        # least two sustained A/C/G-like components before confirming genotype I.
        if assigned_count == 0 or ambiguous_count > len(assignments) * 0.60:
            return (
                "Indeterminate",
                "Whole-genome reference comparison nominated genotype I, but too many parental-window assignments were ambiguous to confirm the expected A/C/G-like mosaic.",
                "",
            )

        i_components = [g for g in strong if g in GENOTYPE_I_PARENTAL_SET]
        i_text = "/".join(i_components)
        if len(i_components) >= 2:
            return (
                "Genotype-I-compatible mosaic",
                f"Whole-genome reference comparison nominated genotype I and sustained parental-like window components ({i_text}) independently support the expected A/C/G-like mosaic.",
                i_text,
            )

        dominant = counts.most_common(1)[0][0] if counts else ""
        dominant_text = f" The parental-window profile was dominated by genotype {dominant}." if dominant else ""
        return (
            "Genotype I not confirmed",
            "Whole-genome reference comparison nominated genotype I, but the independent parental-window screen did not resolve at least two sustained A/C/G-like components." + dominant_text,
            i_text or dominant,
        )

    if assigned_count == 0 or ambiguous_count > len(assignments) * 0.60:
        return (
            "Indeterminate",
            "Too many windows were ambiguous for a reliable recombination screen.",
            "",
        )

    # For conventional A-H genotypes, a recombinant flag requires a sustained
    # *alternative* component. The overall genotype itself does not count as
    # evidence of mosaicism. By default, the alternative must occupy >=20% of
    # informative windows and persist for >=4 consecutive overlapping windows.
    alternatives = [g for g in strong if g != whole_genotype]
    if alternatives:
        components = [whole_genotype] if whole_genotype else []
        components.extend(g for g in alternatives if g not in components)
        components_text = "/".join(components)
        details = ", ".join(
            f"{g}: {counts[g]}/{assigned_count} windows, longest run {longest_runs.get(g, 0)}"
            for g in alternatives
        )
        return (
            "Possible intergenotypic recombination",
            f"A sustained parental-like component discordant with the whole-genome genotype was detected ({details}). Formal recombination/phylogenetic confirmation is recommended.",
            components_text,
        )

    dominant = counts.most_common(1)[0][0] if counts else ""
    return (
        "No significant intergenotypic mosaic detected",
        f"Windowed screening was dominated by genotype {dominant or whole_genotype or 'one genotype'} and no alternative genotype met the minimum sustained-component criteria.",
        dominant or whole_genotype,
    )


def _downgrade_unconfirmed_genotype_i(result, whole, reason):
    """Do not report genotype I unless the independent mosaic screen confirms it."""
    if str(whole.get("genotype", "") or "").upper() != "I":
        return
    support = str(result.get("hbv_curated_genotype_support", "") or "").strip()
    note = (
        "Genotype I was the closest whole-genome reference candidate but was not "
        "independently confirmed by the parental-window A/C/G mosaic screen"
    )
    if reason:
        note += f": {reason}"
    note += "."
    result.update({
        "hbv_curated_genotype": "",
        "hbv_curated_subgenotype": "",
        "hbv_curated_genotype_display": "Indeterminate (I candidate)",
        "hbv_curated_genotype_confidence": "Low",
        "hbv_curated_genotype_support": (support + " " + note).strip(),
    })


def analyse_hbv_genotype_and_recombination(fasta_path, threads=1, timeout=900):
    """Return curated whole-genome HBV genotype + recombination-screen fields."""
    result = {
        "hbv_genotyping_status": "Not run",
        "hbv_genotyping_method": "Whole-genome HBV reference BLAST + conditional sliding-window parental-genotype screen",
        "hbv_genotyping_database": "",
        "hbv_curated_genotype": "",
        "hbv_curated_subgenotype": "",
        "hbv_curated_genotype_display": "",
        "hbv_curated_genotype_confidence": "",
        "hbv_curated_genotype_support": "",
        "hbv_curated_reference_accession": "",
        "hbv_curated_reference_identity": "",
        "hbv_curated_reference_query_coverage": "",
        "hbv_genotype_callable_bases": 0,
        "hbv_genotype_callable_percent": 0.0,
        "hbv_genotype_largest_n_gap": 0,
        "hbv_recombination_status": "Not run",
        "hbv_recombination_summary": "",
        "hbv_recombination_components": "",
        "hbv_recombination_window_size": "",
        "hbv_recombination_step_size": "",
        "hbv_recombination_windows": [],
    }

    db_prefix = find_hbv_genotype_database()
    if not db_prefix:
        result.update({
            "hbv_genotyping_status": "Unavailable",
            "hbv_recombination_status": "Not run",
            "hbv_curated_genotype_support": (
                f"HBV genotype reference database not found. Set {HBV_GENOTYPE_DB_ENV} or build the package-local database."
            ),
        })
        return result

    blastn = shutil.which("blastn")
    if not blastn:
        result.update({
            "hbv_genotyping_status": "Unavailable",
            "hbv_curated_genotype_support": "blastn is not available on PATH.",
        })
        return result

    result["hbv_genotyping_database"] = db_prefix
    try:
        proc = _run_blast(blastn, fasta_path, db_prefix, threads, timeout, task="megablast", max_targets=100)
    except subprocess.TimeoutExpired:
        result.update({"hbv_genotyping_status": "Timeout", "hbv_curated_genotype_support": "HBV whole-genome genotype BLAST timed out."})
        return result
    except Exception as exc:
        result.update({"hbv_genotyping_status": "Error", "hbv_curated_genotype_support": str(exc)})
        return result

    if proc.returncode != 0:
        result.update({
            "hbv_genotyping_status": "Error",
            "hbv_curated_genotype_support": (proc.stderr or "HBV genotype BLAST failed").strip(),
        })
        return result

    sequence = _read_single_fasta_sequence(fasta_path)
    sequence_quality = _sequence_quality(sequence)
    whole_hits = _parse_output(proc.stdout or "")
    whole = _whole_genome_call(whole_hits, sequence_quality)
    best = whole.get("best_hit") or {}
    result.update({
        "hbv_genotyping_status": "OK",
        "hbv_curated_genotype": whole.get("genotype", ""),
        "hbv_curated_subgenotype": whole.get("subgenotype", ""),
        "hbv_curated_genotype_display": whole.get("display", ""),
        "hbv_curated_genotype_confidence": whole.get("confidence", ""),
        "hbv_curated_genotype_support": whole.get("support", ""),
        "hbv_curated_reference_accession": best.get("source_accession", "") or best.get("subject_id", ""),
        "hbv_curated_reference_identity": round(float(best.get("pident", 0) or 0), 3) if best else "",
        "hbv_curated_reference_query_coverage": round(float(best.get("qcovs", 0) or 0), 1) if best else "",
        "hbv_genotype_callable_bases": sequence_quality.get("callable_bases", 0),
        "hbv_genotype_callable_percent": round(float(sequence_quality.get("callable_percent", 0.0) or 0.0), 1),
        "hbv_genotype_largest_n_gap": sequence_quality.get("largest_n_gap", 0),
    })

    # Recombination is a whole-genome pattern. Do not infer presence or absence
    # of mosaicism when too much of the genome is missing or a large contiguous
    # gap could hide a breakpoint/component. Genotyping can still be reportable.
    callable_percent = float(sequence_quality.get("callable_percent", 0.0) or 0.0)
    largest_n_gap = int(sequence_quality.get("largest_n_gap", 0) or 0)
    if callable_percent < 80.0 or largest_n_gap > 500:
        result.update({
            "hbv_recombination_status": "Not evaluable",
            "hbv_recombination_summary": (
                f"Rekombinationsbedömning ej möjlig: endast {callable_percent:.1f}% av konsensus är analyserbar "
                f"och största sammanhängande N-gap är {largest_n_gap} nt. Genotypning baseras endast på tillgänglig sekvens."
            ),
            "hbv_recombination_components": "",
            "hbv_recombination_windows": [],
        })
        _downgrade_unconfirmed_genotype_i(result, whole, result["hbv_recombination_summary"])
        return result

    indicated, trigger_reason = _recombination_indicated(whole)
    if not indicated:
        result.update({
            "hbv_recombination_status": "Not indicated",
            "hbv_recombination_summary": (
                "Recombination screening was not run: " + trigger_reason
            ),
            "hbv_recombination_components": whole.get("genotype", ""),
            "hbv_recombination_windows": [],
        })
        return result

    window = _env_int(HBV_RECOMB_WINDOW_ENV, DEFAULT_WINDOW, minimum=250)
    step = _env_int(HBV_RECOMB_STEP_ENV, DEFAULT_STEP, minimum=50)
    min_identity = _env_float(HBV_RECOMB_MIN_IDENTITY_ENV, DEFAULT_MIN_IDENTITY)
    min_qcov = _env_float(HBV_RECOMB_MIN_QCOV_ENV, DEFAULT_MIN_QCOV)
    min_margin = _env_float(HBV_RECOMB_MIN_MARGIN_ENV, DEFAULT_MIN_IDENTITY_MARGIN)
    windows = _make_windows(sequence, window, step)
    result["hbv_recombination_window_size"] = window
    result["hbv_recombination_step_size"] = step

    if not windows:
        result.update({
            "hbv_recombination_status": "Indeterminate",
            "hbv_recombination_summary": "Consensus is too short for the configured recombination window screen.",
        })
        _downgrade_unconfirmed_genotype_i(result, whole, result["hbv_recombination_summary"])
        return result

    with tempfile.TemporaryDirectory(prefix="aurora_hbv_recomb_") as tmpdir:
        query_path = os.path.join(tmpdir, "windows.fasta")
        with open(query_path, "w", encoding="utf-8") as handle:
            for item in windows:
                if item["n_fraction"] <= 0.20:
                    handle.write(f">{item['id']}\n{item['sequence']}\n")
        try:
            proc = _run_blast(blastn, query_path, db_prefix, threads, timeout, task="blastn", max_targets=100)
        except subprocess.TimeoutExpired:
            result.update({
                "hbv_recombination_status": "Indeterminate",
                "hbv_recombination_summary": "Sliding-window HBV recombination BLAST timed out.",
            })
            _downgrade_unconfirmed_genotype_i(result, whole, result["hbv_recombination_summary"])
            return result
        except Exception as exc:
            result.update({
                "hbv_recombination_status": "Indeterminate",
                "hbv_recombination_summary": f"Sliding-window screen failed: {exc}",
            })
            _downgrade_unconfirmed_genotype_i(result, whole, result["hbv_recombination_summary"])
            return result

    if proc.returncode != 0:
        result.update({
            "hbv_recombination_status": "Indeterminate",
            "hbv_recombination_summary": (proc.stderr or "Sliding-window HBV BLAST failed").strip(),
        })
        _downgrade_unconfirmed_genotype_i(result, whole, result["hbv_recombination_summary"])
        return result

    by_query = defaultdict(list)
    for hit in _parse_output(proc.stdout or ""):
        by_query[hit.get("qseqid", "")].append(hit)

    assignments = []
    for window_item in windows:
        if window_item["n_fraction"] > 0.20:
            assignments.append({
                **window_item,
                "genotype": "Ambiguous",
                "identity": "",
                "qcov": "",
                "margin": "",
                "accession": "",
            })
        else:
            assignments.append(_assign_window(
                window_item,
                by_query.get(window_item["id"], []),
                min_identity,
                min_qcov,
                min_margin,
            ))

    recomb_status, recomb_summary, components = _classify_recombination(
        assignments, whole.get("genotype", "")
    )
    result.update({
        "hbv_recombination_status": recomb_status,
        "hbv_recombination_summary": recomb_summary,
        "hbv_recombination_components": components,
        "hbv_recombination_windows": [
            {k: v for k, v in item.items() if k != "sequence"}
            for item in assignments
        ],
    })
    if str(whole.get("genotype", "") or "").upper() == "I" and recomb_status != "Genotype-I-compatible mosaic":
        _downgrade_unconfirmed_genotype_i(result, whole, recomb_summary)
    return result
