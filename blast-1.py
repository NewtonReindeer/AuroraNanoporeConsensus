"""Local NCBI BLAST integration and conservative virus genotype inference."""
from __future__ import annotations

from .cancellation import run_cancellable
from .hbv_genotyping import analyse_hbv_genotype_and_recombination

import csv
import os
import re
import shutil
import subprocess
from collections import Counter
from datetime import datetime
from io import StringIO


BLAST_DB_ENV = "AURORA_BLAST_DB"
BLAST_SMALL_DB_ROOT_ENV = "AURORA_BLAST_SMALL_DB_ROOT"
BLAST_THREADS_ENV = "AURORA_BLAST_THREADS"
BLAST_TIMEOUT_ENV = "AURORA_BLAST_TIMEOUT_SECONDS"
DEFAULT_MAX_TARGET_SEQS = 25

# Dedicated Aurora BLAST databases, relative to the small-database root.
# These are preferred because they are dramatically faster than searching the
# complete nt_viruses database with a taxid filter.
SMALL_DATABASE_RELATIVE_PATHS = {
    "HBV": os.path.join("hbv", "hbv"),
    "HCV": os.path.join("hcv", "hcv"),
    "ADENOVIRUS": os.path.join("adenovirus", "adenovirus"),
    "CMV": os.path.join("cmv", "cmv"),
}

# Taxonomy restrictions are used only as a fallback when the matching small
# database is unavailable. "Other" intentionally has no restriction.
DEFAULT_VIRUS_TAXIDS = {
    "HBV": "10407",
    "HCV": "3052230",
    "ADENOVIRUS": "129875,108098,129951,130310,130308,130309,536079",
    "CMV": "3050295",
}


# Curated accession overrides for HBV records whose GenBank title does not
# reliably expose the accepted genotype.  Keep this deliberately small and
# evidence-based; it is consulted before any title-based inference.
HBV_ACCESSION_GENOTYPE = {
    "HQ231884": ("I", ""),
}


def _normalise_accession(accession):
    """Return an accession without version suffix, e.g. HQ231884.1 -> HQ231884."""
    return str(accession or "").strip().split(".", 1)[0].upper()

BLAST_TAXID_ENV = {
    "HBV": "AURORA_BLAST_HBV_TAXIDS",
    "HCV": "AURORA_BLAST_HCV_TAXIDS",
    "ADENOVIRUS": "AURORA_BLAST_ADENOVIRUS_TAXIDS",
    "CMV": "AURORA_BLAST_CMV_TAXIDS",
}


def _package_dir():
    return os.path.dirname(os.path.abspath(__file__))


def _blast_db_candidates():
    """Yield candidate prefixes for the full local nt_viruses database."""
    candidates = []
    explicit = os.environ.get(BLAST_DB_ENV, "").strip()
    if explicit:
        candidates.append(explicit)

    # Development/package-local layout:
    # aurora_refactored/blastdb/nt_viruses
    candidates.append(os.path.join(_package_dir(), "blastdb", "nt_viruses"))

    candidates.extend([
        "/mnt/NewtVault/AdenoPipeline/Database/nt_viruses",
        "/media/NewtVault/AdenoPipeline/Database/nt_viruses",
        "/data/AdenoPipeline/Database/nt_viruses",
        "/opt/aurora/Database/nt_viruses",
        "/Volumes/NewtVault/AdenoPipeline/Database/nt_viruses",
        "/Volumes/Newt-Vault/AdenoPipeline/Database/nt_viruses",
    ])
    seen = set()
    for path in candidates:
        if path and path not in seen:
            seen.add(path)
            yield path


def _small_db_root_candidates(master_db_prefix=""):
    """Yield candidate roots containing hbv/, hcv/, adenovirus/, and cmv/."""
    candidates = []
    explicit = os.environ.get(BLAST_SMALL_DB_ROOT_ENV, "").strip()
    if explicit:
        candidates.append(explicit)

    # Preferred development/package-local layout:
    # aurora_refactored/blastdb/aurora_small_db/...
    candidates.append(os.path.join(_package_dir(), "blastdb", "aurora_small_db"))

    # If AURORA_BLAST_DB points at .../blastdb/nt_viruses, automatically look
    # beside it for .../blastdb/aurora_small_db.
    if master_db_prefix:
        candidates.append(os.path.join(os.path.dirname(master_db_prefix), "aurora_small_db"))

    candidates.extend([
        "/mnt/NewtVault/AdenoPipeline/Database/aurora_small_db",
        "/media/NewtVault/AdenoPipeline/Database/aurora_small_db",
        "/data/AdenoPipeline/Database/aurora_small_db",
        "/opt/aurora/Database/aurora_small_db",
        "/Volumes/NewtVault/AdenoPipeline/Database/aurora_small_db",
        "/Volumes/Newt-Vault/AdenoPipeline/Database/aurora_small_db",
    ])

    seen = set()
    for path in candidates:
        if path and path not in seen:
            seen.add(path)
            yield path


def _db_prefix_exists(prefix):
    if not prefix:
        return False
    suffixes = (".nal", ".ndb", ".nin", ".00.nin", ".nsq", ".00.nsq")
    return any(os.path.exists(prefix + suffix) for suffix in suffixes)


def _validate_blast_database(candidate):
    if not _db_prefix_exists(candidate):
        return False
    blastdbcmd = shutil.which("blastdbcmd")
    if not blastdbcmd:
        return True
    try:
        proc = run_cancellable(
            [blastdbcmd, "-db", candidate, "-info"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=30,
            check=False,
        )
        return proc.returncode == 0
    except Exception:
        return False


def find_blast_database():
    """Return a usable full nt_viruses database prefix, or an empty string."""
    for candidate in _blast_db_candidates():
        if _validate_blast_database(candidate):
            return candidate
    return ""


def find_small_blast_database(virus_type, master_db_prefix=""):
    """Return the dedicated small database prefix for a known virus, if present."""
    virus = _normalise_virus_type(virus_type)
    relative = SMALL_DATABASE_RELATIVE_PATHS.get(virus)
    if not relative:
        return ""

    for root in _small_db_root_candidates(master_db_prefix):
        candidate = os.path.join(root, relative)
        if _validate_blast_database(candidate):
            return candidate
    return ""


def blast_database_info(db_prefix):
    result = {"blast_database": db_prefix or "", "blast_database_date": "", "blast_database_sequences": ""}
    blastdbcmd = shutil.which("blastdbcmd")
    if not db_prefix or not blastdbcmd:
        return result
    try:
        proc = run_cancellable(
            [blastdbcmd, "-db", db_prefix, "-info"],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if proc.returncode != 0:
            return result
        text = proc.stdout or ""
        date_match = re.search(r"^Date:\s*(.+?)\s+Longest sequence:", text, re.M)
        seq_match = re.search(r"([\d,]+)\s+sequences;", text)
        if date_match:
            result["blast_database_date"] = date_match.group(1).strip()
        if seq_match:
            result["blast_database_sequences"] = seq_match.group(1).replace(",", "")
    except Exception:
        pass
    return result


def empty_blast_result(status="Not run", message=""):
    return {
        "blast_status": status,
        "blast_message": message,
        "blast_database": "",
        "blast_database_date": "",
        "blast_database_sequences": "",
        "blast_search_scope": "",
        "blast_taxid_filter": "",
        "blast_top_accession": "",
        "blast_top_organism": "",
        "blast_top_title": "",
        "blast_top_identity": "",
        "blast_top_query_coverage": "",
        "blast_top_evalue": "",
        "blast_top_bitscore": "",
        "blast_genotype": "",
        "blast_subtype": "",
        "blast_genotype_display": "",
        "blast_genotype_confidence": "",
        "blast_genotype_support": "",
        "blast_genotype_match_accession": "",
        "blast_genotype_match_organism": "",
        "blast_genotype_match_title": "",
        "blast_genotype_match_identity": "",
        "blast_genotype_match_query_coverage": "",
        "blast_genotype_match_evalue": "",
        "blast_genotype_match_bitscore": "",
        "blast_hits": [],
        "blast_run_time": "",
        "hbv_genotyping_status": "",
        "hbv_genotyping_method": "",
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
        "hbv_recombination_status": "",
        "hbv_recombination_summary": "",
        "hbv_recombination_components": "",
        "hbv_recombination_window_size": "",
        "hbv_recombination_step_size": "",
        "hbv_recombination_windows": [],
    }


def _normalise_virus_type(virus_type):
    text = str(virus_type or "").strip().lower()
    if "hepatitis b" in text or text in {"hbv", "hepatit b"}:
        return "HBV"
    if "hepatitis c" in text or text in {"hcv", "hepatit c"}:
        return "HCV"
    if "adeno" in text:
        return "ADENOVIRUS"
    if "cytomegalo" in text or text in {"cmv", "hcmv", "human cmv"}:
        return "CMV"
    if text == "other":
        return "OTHER"
    return text.upper() or "VIRUS"




def _taxids_for_virus(virus_type):
    """Return the BLAST taxid restriction for a selected Aurora virus pipeline.

    HBV, HCV, Adenovirus, and CMV use a restricted full-database fallback.
    The Other pipeline returns an empty string, which means the complete local
    nt_viruses database is searched.
    Environment variables can override the defaults without changing source code.
    """
    virus = _normalise_virus_type(virus_type)
    default = DEFAULT_VIRUS_TAXIDS.get(virus, "")
    env_name = BLAST_TAXID_ENV.get(virus)
    if env_name:
        return os.environ.get(env_name, default).strip()
    return ""


def _search_scope_for_virus(virus_type, taxids, dedicated_db=False):
    virus = _normalise_virus_type(virus_type)
    display_names = {
        "HBV": "HBV",
        "HCV": "HCV",
        "ADENOVIRUS": "human adenovirus",
        "CMV": "human CMV",
    }
    label = display_names.get(virus, virus)
    if dedicated_db:
        return f"Dedicated {label} database"
    if not taxids:
        return "All viruses (unrestricted nt_viruses)"
    return f"{label}-restricted nt_viruses fallback"

def _extract_hbv_genotype(title):
    """Extract HBV genotype information from a BLAST subject title.

    Explicit ``genotype D`` / ``subgenotype D2`` annotations are preferred.
    Some GenBank HBV records instead encode the genotype in an isolate or
    strain name, for example ``isolate D_NZL_AN10_1984``.  That second form is
    accepted as weaker supporting evidence, but only when the genotype letter
    is immediately followed by ``_`` or ``-`` after the isolate/strain label.

    Returns ``(genotype, subtype, evidence_type)`` where evidence_type is
    ``explicit`` or ``isolate_prefix``.
    """
    text = str(title or "")

    sub = re.search(r"\bsubgenotype\s*[:=]?\s*([A-J](?:\d+)?)\b", text, re.I)
    geno = re.search(r"\bgenotype\s*[:=]?\s*([A-J])(?:\b|\d)", text, re.I)
    subtype = sub.group(1).upper() if sub else ""
    genotype = geno.group(1).upper() if geno else (subtype[:1] if subtype else "")
    if genotype:
        return genotype, subtype, "explicit"

    # Conservative fallback for titles such as:
    #   Hepatitis B virus isolate D_NZL_AN10_1984, complete genome
    # Do not interpret an arbitrary standalone A-J elsewhere in the title.
    isolate = re.search(
        r"\b(?:isolate|strain)\s*[:=]?\s*([A-J])[_-][A-Za-z0-9]",
        text,
        re.I,
    )
    if isolate:
        return isolate.group(1).upper(), "", "isolate_prefix"

    return "", "", ""


def _extract_hcv_genotype(title):
    text = str(title or "")
    patterns = [
        r"\bsubtype\s*[:=]?\s*([1-8][a-z]?)\b",
        r"\bgenotype\s*[:=]?\s*([1-8](?:[a-z])?)\b",
        r"\bgenotype\s*([1-8])\s*subtype\s*([a-z])\b",
    ]
    for idx, pattern in enumerate(patterns):
        match = re.search(pattern, text, re.I)
        if not match:
            continue
        if idx == 2:
            subtype = (match.group(1) + match.group(2)).lower()
        else:
            subtype = match.group(1).lower()
        genotype = re.match(r"[1-8]", subtype).group(0)
        return genotype, subtype
    return "", ""


def _infer_genotype(virus_type, hits):
    """Infer genotype from high-quality BLAST hits.

    HBV deliberately does *not* use majority voting across GenBank titles.
    First, an exact accession override is used for a small set of curated
    records whose title is not reliable for genotype inference.  Otherwise,
    the highest-ranking hit that itself carries genotype evidence is used.

    HCV keeps Aurora's existing multi-hit consensus logic unchanged.
    """
    virus = _normalise_virus_type(virus_type)
    if not hits:
        return "", "", "", "", ""

    best_identity = max(float(h.get("pident", 0) or 0) for h in hits)

    if virus == "HBV":
        # 1) Exact curated accession takes precedence. This fixes cases such as
        # HQ231884.1 (genotype I), whose BLAST title may not contain a usable
        # genotype annotation.
        for hit in hits[:25]:
            identity = float(hit.get("pident", 0) or 0)
            qcov = float(hit.get("qcovs", 0) or 0)
            if qcov < 70.0 or identity < max(80.0, best_identity - 2.0):
                continue
            accession_key = _normalise_accession(hit.get("accession", ""))
            mapped = HBV_ACCESSION_GENOTYPE.get(accession_key)
            if mapped:
                genotype, subtype = mapped
                # Genotype I is a recognized recombinant lineage and must not be
                # established from a single legacy BLAST accession override.
                # The HBV-specific curated classifier below performs the required
                # independent parental-window confirmation when available.
                if genotype == "I":
                    continue
                display = genotype if not subtype or subtype == genotype else f"{genotype} ({subtype})"
                support = (
                    f"HBV genotype from reference accession {hit.get('accession', '')}; "
                    f"identity {identity:.3f}%, query coverage {qcov:.1f}%"
                )
                confidence = "High" if qcov >= 90.0 and identity >= 90.0 else "Moderate"
                return genotype, subtype, display, confidence, support

        # 2) No accession override: use the strongest individual hit that
        # actually carries genotype evidence. Do not count titles.
        candidates = []
        for hit in hits[:25]:
            identity = float(hit.get("pident", 0) or 0)
            qcov = float(hit.get("qcovs", 0) or 0)
            if qcov < 70.0 or identity < max(80.0, best_identity - 2.0):
                continue
            genotype, subtype, evidence_type = _extract_hbv_genotype(hit.get("stitle", ""))
            if not genotype:
                continue
            evidence_rank = 1 if evidence_type == "explicit" else 0
            rank = (
                evidence_rank,
                float(hit.get("bitscore", 0) or 0),
                qcov,
                identity,
            )
            candidates.append((rank, hit, genotype, subtype, evidence_type))

        if not candidates:
            return "", "", "", "Low", "No reliable HBV genotype annotation found among high-quality BLAST hits"

        candidates.sort(key=lambda item: item[0], reverse=True)
        _, best_hit, genotype, subtype, evidence_type = candidates[0]
        if genotype == "I":
            return (
                "", "", "Indeterminate (I candidate)", "Low",
                "Legacy BLAST nominated genotype I, but genotype I requires independent confirmation by the HBV-specific parental-window screen."
            )
        display = genotype if not subtype or subtype == genotype else f"{genotype} ({subtype})"
        support = (
            f"HBV genotype from strongest {evidence_type} BLAST hit "
            f"{best_hit.get('accession', '')}; identity "
            f"{float(best_hit.get('pident', 0) or 0):.3f}%, query coverage "
            f"{float(best_hit.get('qcovs', 0) or 0):.1f}%"
        )
        confidence = "High" if evidence_type == "explicit" and float(best_hit.get("qcovs", 0) or 0) >= 90.0 else "Moderate"
        return genotype, subtype, display, confidence, support

    # HCV: preserve the previous multi-hit consensus behaviour.
    parsed = []
    for hit in hits[:25]:
        identity = float(hit.get("pident", 0) or 0)
        qcov = float(hit.get("qcovs", 0) or 0)
        if qcov < 70.0 or identity < max(80.0, best_identity - 2.0):
            continue
        title = hit.get("stitle", "")
        if virus == "HCV":
            genotype, subtype = _extract_hcv_genotype(title)
        else:
            genotype, subtype = "", ""
        if genotype:
            parsed.append((genotype, subtype, identity, qcov))

    if not parsed:
        return "", "", "", "", ""

    total = len(parsed)
    genotype_counts = Counter(item[0] for item in parsed)
    genotype, genotype_count = genotype_counts.most_common(1)[0]
    genotype_fraction = genotype_count / total if total else 0.0
    matching_genotype = [item for item in parsed if item[0] == genotype]

    subtype = ""
    subtype_count = 0
    subtype_values = [item[1] for item in matching_genotype if item[1]]
    if subtype_values:
        subtype, subtype_count = Counter(subtype_values).most_common(1)[0]

    if virus == "HCV":
        if subtype and subtype_count / max(len(matching_genotype), 1) >= 0.60:
            display = subtype
            support = f"{subtype_count}/{total} parseable top hits support subtype {subtype}"
        else:
            display = genotype
            support = f"{genotype_count}/{total} parseable top hits support genotype {genotype}"

        if genotype_count >= 3 and genotype_fraction >= 0.80:
            confidence = "High"
        elif genotype_count >= 2 and genotype_fraction >= 0.60:
            confidence = "Moderate"
        else:
            confidence = "Low"
    else:
        display = subtype or genotype
        support = f"{genotype_count}/{total} parseable top hits"
        confidence = "Low"

    return genotype, subtype, display, confidence, support

def _best_genotype_match(virus_type, hits, genotype, subtype=""):
    """Return the best BLAST hit that actually carries genotype evidence.

    The overall highest-scoring BLAST hit can be an untyped isolate, a gene, or
    a partial record.  For report display we instead prefer a hit whose title
    supports the inferred genotype.  Complete-genome records are preferred over
    partial/gene records, then query coverage, identity, and bit score are used
    to choose the closest match.  The original overall top hit remains stored in
    ``blast_top_*`` for traceability.
    """
    virus = _normalise_virus_type(virus_type)
    if not hits or not genotype:
        return None

    candidates = []
    for hit in hits:
        title = str(hit.get("stitle", "") or "")
        evidence_type = ""
        hit_genotype = ""
        hit_subtype = ""

        if virus == "HBV":
            hit_genotype, hit_subtype, evidence_type = _extract_hbv_genotype(title)
        elif virus == "HCV":
            hit_genotype, hit_subtype = _extract_hcv_genotype(title)
            evidence_type = "explicit" if hit_genotype else ""
        else:
            continue

        if hit_genotype != genotype:
            continue

        # If a subtype was confidently inferred, prefer the same subtype when
        # available, but do not discard genotype-only records.
        subtype_match = bool(subtype and hit_subtype and hit_subtype.lower() == str(subtype).lower())
        explicit = evidence_type == "explicit"
        complete_genome = bool(re.search(r"\bcomplete genome\b", title, re.I))
        partial_penalty = bool(re.search(
            r"\b(?:partial|fragment|gene|protein|precore|pre-core|core protein|surface antigen)\b",
            title, re.I,
        ))

        score = (
            1 if subtype_match else 0,
            1 if explicit else 0,
            1 if complete_genome else 0,
            0 if partial_penalty else 1,
            float(hit.get("qcovs", 0) or 0),
            float(hit.get("pident", 0) or 0),
            float(hit.get("bitscore", 0) or 0),
        )
        candidates.append((score, hit))

    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def _organism_from_title(title):
    text = str(title or "").strip()
    if not text:
        return ""
    for known in ("Hepatitis B virus", "Hepatitis C virus", "Human cytomegalovirus"):
        if known.lower() in text.lower():
            return known
    adeno = re.search(r"\b(?:Human )?(?:mastadenovirus|adenovirus)\s+[A-Za-z0-9-]+", text, re.I)
    if adeno:
        return adeno.group(0)
    match = re.match(r"^(.+?)(?:,|\s+(?:isolate|strain|genotype|subtype|complete genome|segment)\b|$)", text, re.I)
    return (match.group(1) if match else text[:100]).strip()

def _parse_blast_output(text):
    hits = []
    if not text.strip():
        return hits
    reader = csv.reader(StringIO(text), delimiter="\t")
    for fields in reader:
        if len(fields) < 9:
            continue
        try:
            hit = {
                "qseqid": fields[0],
                "accession": fields[1],
                "pident": float(fields[2]),
                "qcovs": float(fields[3]),
                "length": int(float(fields[4])),
                "evalue": fields[5],
                "bitscore": float(fields[6]),
                "taxids": fields[7],
                "stitle": fields[8],
                "organism": _organism_from_title(fields[8]),
            }
        except (ValueError, IndexError):
            continue
        hits.append(hit)
    return hits


def _run_local_blast(
    blastn, fasta_path, db_prefix, threads, timeout,
    max_target_seqs=DEFAULT_MAX_TARGET_SEQS, taxids="",
):
    """Run one local megablast search and return (process, elapsed_seconds)."""
    outfmt = "6 qseqid saccver pident qcovs length evalue bitscore staxids stitle"
    cmd = [
        blastn,
        "-task", "megablast",
        "-query", fasta_path,
        "-db", db_prefix,
        "-max_target_seqs", str(max_target_seqs),
        "-evalue", "1e-20",
        "-num_threads", str(threads),
        "-outfmt", outfmt,
    ]
    if taxids:
        cmd.extend(["-taxids", taxids])

    started = datetime.now()
    proc = run_cancellable(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    return proc, (datetime.now() - started).total_seconds()


def blast_consensus(fasta_path, virus_type, max_target_seqs=DEFAULT_MAX_TARGET_SEQS):
    """BLAST a Medaka consensus and return fields used by Aurora reports.

    Known-virus pipelines prefer the small derived database for speed. If that
    database is missing, Aurora falls back to the full nt_viruses database with
    the appropriate taxid restriction. The Other pipeline always searches the
    complete nt_viruses database without a taxonomy restriction.
    """
    result = empty_blast_result()
    if not fasta_path or not os.path.isfile(fasta_path) or os.path.getsize(fasta_path) == 0:
        result.update({"blast_status": "Skipped", "blast_message": "Consensus FASTA missing or empty."})
        return result

    blastn = shutil.which("blastn")
    if not blastn:
        print("[BLAST ERROR] blastn is not available on PATH.")
        result.update({"blast_status": "Unavailable", "blast_message": "blastn is not available on PATH."})
        return result

    virus = _normalise_virus_type(virus_type)
    master_db = find_blast_database()
    small_db = find_small_blast_database(virus_type, master_db)

    # Known viruses: prefer the dedicated small database.
    # Other: always use the full nt_viruses database.
    if small_db and virus in SMALL_DATABASE_RELATIVE_PATHS:
        db_prefix = small_db
        taxids = ""
        using_dedicated_db = True
    else:
        db_prefix = master_db
        taxids = _taxids_for_virus(virus_type) if virus != "OTHER" else ""
        using_dedicated_db = False

    if not db_prefix:
        print(f"[BLAST ERROR] No local BLAST database for {virus_type}. "
              f"Check {BLAST_DB_ENV} and {BLAST_SMALL_DB_ROOT_ENV}.")
        result.update({
            "blast_status": "Unavailable",
            "blast_message": (
                "No usable local BLAST database was found. "
                f"Set {BLAST_DB_ENV} to the nt_viruses database prefix. "
                f"Optionally set {BLAST_SMALL_DB_ROOT_ENV} to the aurora_small_db directory."
            ),
        })
        return result

    search_scope = _search_scope_for_virus(virus_type, taxids, dedicated_db=using_dedicated_db)
    result["blast_search_scope"] = search_scope
    result["blast_taxid_filter"] = taxids

    db_info = blast_database_info(db_prefix)
    threads_default = min(max(os.cpu_count() or 1, 1), 8)
    try:
        threads = max(1, int(os.environ.get(BLAST_THREADS_ENV, threads_default)))
    except ValueError:
        threads = threads_default
    try:
        timeout = max(60, int(os.environ.get(BLAST_TIMEOUT_ENV, "900")))
    except ValueError:
        timeout = 900

    if using_dedicated_db:
        print(f"BLAST search scope: {search_scope}")
        print(f"BLAST database: {db_prefix}")
    elif taxids:
        print(f"BLAST search scope: {search_scope} (taxid filter: {taxids})")
        print(f"BLAST database: {db_prefix}")
        print("[BLAST WARNING] Dedicated small database not found; using taxid-filtered nt_viruses fallback.")
    else:
        print(f"BLAST search scope: {search_scope}")
        print(f"BLAST database: {db_prefix}")

    started = datetime.now()
    result["blast_run_time"] = started.isoformat(timespec="seconds")
    try:
        proc, elapsed = _run_local_blast(
            blastn, fasta_path, db_prefix, threads, timeout,
            max_target_seqs=max_target_seqs, taxids=taxids,
        )
    except subprocess.TimeoutExpired:
        result.update(db_info)
        result.update({"blast_status": "Timeout", "blast_message": f"BLAST exceeded {timeout} seconds."})
        return result
    except Exception as exc:
        result.update(db_info)
        result.update({"blast_status": "Error", "blast_message": str(exc)})
        return result

    print(f"BLAST completed in {elapsed:.2f} seconds.")

    result.update(db_info)
    if proc.returncode != 0:
        result.update({"blast_status": "Error", "blast_message": (proc.stderr or "BLAST failed").strip()})
        return result

    hits = _parse_blast_output(proc.stdout or "")
    result["blast_hits"] = hits
    if not hits:
        scope_note = f" within taxid filter {taxids}" if taxids else ""
        result.update({
            "blast_status": "No hits",
            "blast_message": f"No local BLAST hit was returned{scope_note}.",
        })
        return result

    top = hits[0]
    genotype, subtype, display, confidence, support = _infer_genotype(virus_type, hits)
    genotype_match = _best_genotype_match(virus_type, hits, genotype, subtype)

    result.update({
        "blast_status": "OK",
        "blast_message": "",
        "blast_top_accession": top.get("accession", ""),
        "blast_top_organism": top.get("organism", ""),
        "blast_top_title": top.get("stitle", ""),
        "blast_top_identity": round(float(top.get("pident", 0) or 0), 3),
        "blast_top_query_coverage": round(float(top.get("qcovs", 0) or 0), 3),
        "blast_top_evalue": top.get("evalue", ""),
        "blast_top_bitscore": round(float(top.get("bitscore", 0) or 0), 2),
        "blast_genotype": genotype,
        "blast_subtype": subtype,
        "blast_genotype_display": display,
        "blast_genotype_confidence": confidence,
        "blast_genotype_support": support,
        "blast_genotype_match_accession": (genotype_match or {}).get("accession", ""),
        "blast_genotype_match_organism": (genotype_match or {}).get("organism", ""),
        "blast_genotype_match_title": (genotype_match or {}).get("stitle", ""),
        "blast_genotype_match_identity": round(float((genotype_match or {}).get("pident", 0) or 0), 3) if genotype_match else "",
        "blast_genotype_match_query_coverage": round(float((genotype_match or {}).get("qcovs", 0) or 0), 3) if genotype_match else "",
        "blast_genotype_match_evalue": (genotype_match or {}).get("evalue", ""),
        "blast_genotype_match_bitscore": round(float((genotype_match or {}).get("bitscore", 0) or 0), 2) if genotype_match else "",
    })

    # HBV gets a second, Aurora-controlled whole-genome genotype comparison plus
    # a conservative sliding-window recombination screen. The ordinary local
    # NCBI BLAST above is retained for identification/traceability.
    if virus == "HBV":
        print("Running HBV genotyping (conditional recombination screen)...")
        hbv_result = analyse_hbv_genotype_and_recombination(
            fasta_path, threads=threads, timeout=timeout
        )
        result.update(hbv_result)
        if hbv_result.get("hbv_genotyping_status") == "OK":
            curated_display = hbv_result.get("hbv_curated_genotype_display", "")
            curated_genotype = hbv_result.get("hbv_curated_genotype", "")
            curated_subgenotype = hbv_result.get("hbv_curated_subgenotype", "")
            # An ambiguous curated result deliberately overrides title-based
            # inference rather than forcing a potentially wrong genotype.
            result["blast_genotype"] = curated_genotype
            result["blast_subtype"] = curated_subgenotype
            result["blast_genotype_display"] = curated_display or "Indeterminate"
            result["blast_genotype_confidence"] = hbv_result.get("hbv_curated_genotype_confidence", "Low")
            result["blast_genotype_support"] = hbv_result.get("hbv_curated_genotype_support", "")

            # These fields now describe the curated genotype reference for HBV.
            result["blast_genotype_match_accession"] = hbv_result.get("hbv_curated_reference_accession", "")
            result["blast_genotype_match_organism"] = "Hepatitis B virus"
            result["blast_genotype_match_title"] = "Aurora HBV genotype reference"
            result["blast_genotype_match_identity"] = hbv_result.get("hbv_curated_reference_identity", "")
            result["blast_genotype_match_query_coverage"] = hbv_result.get("hbv_curated_reference_query_coverage", "")
            result["blast_genotype_match_evalue"] = ""
            result["blast_genotype_match_bitscore"] = ""
        elif hbv_result.get("hbv_genotyping_status") in {"Unavailable", "Error", "Timeout"}:
            # Keep the existing title/top-hit inference as a fallback, but make
            # the absence of curated genotyping explicit in the result fields.
            print("[HBV GENOTYPING WARNING] HBV genotype reference screen unavailable; using legacy HBV BLAST inference as fallback.")

    return result
