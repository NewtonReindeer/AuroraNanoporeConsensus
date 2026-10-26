"""Independent CMV-only BCFtools review branch.

Runs on original FASTQ against fixed Merlin coordinates, never modifies the
primary Aurora/Medaka mutation calls. This is not a validated viral caller.
"""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

ACCESSION = "NC_006273.2"


def _resolve_reference(reference_dir):
    explicit = os.environ.get("AURORA_CMV_MERLIN_FASTA", "").strip()
    folder = Path(reference_dir)
    candidates = ([Path(explicit)] if explicit else
                  [folder / "HumanHerpesvirus5GenomeMerlin.fa",
                   folder / "HumanHerpesvirus5GenomeMerlin.fasta",
                   folder / "Merlin.fa", folder / "Merlin.fasta"])
    for candidate in candidates:
        if candidate.is_file():
            with candidate.open() as handle:
                first = handle.readline().strip()
            if first.split()[0] != ">" + ACCESSION:
                raise ValueError(f"Merlin FASTA header must start >{ACCESSION}; found {first!r} in {candidate}")
            return candidate.resolve()
    raise FileNotFoundError("Merlin FASTA not found. Set AURORA_CMV_MERLIN_FASTA to the NC_006273.2 FASTA path")


def _checked(argv, log, **kwargs):
    log.write("COMMAND: " + " ".join(map(str, argv)) + "\n")
    log.flush()
    subprocess.run(list(map(str, argv)), check=True, stderr=log, **kwargs)


def _vcf_to_tsv(bcftools, vcf, tsv, log):
    cmd = [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\n", str(vcf)]
    with tsv.open("w") as out:
        out.write("chrom\tposition\tref\talt\tqual\tfilter\n")
        _checked(cmd, log, stdout=out)
    with tsv.open() as handle:
        return max(0, sum(1 for _ in handle) - 1)


def run(sample, fastq, reference_dir, sample_dir):
    if not sample or Path(sample).name != sample or "/" in sample or "\\" in sample:
        raise ValueError("Invalid sample ID")
    fastq = Path(fastq).resolve()
    if not fastq.is_file():
        raise FileNotFoundError(f"FASTQ not found: {fastq}")
    reference = _resolve_reference(reference_dir)
    tools = {x: shutil.which(x) for x in ("minimap2", "samtools", "bcftools")}

    # Prefer BCFtools from the active Conda environment
    import sys
    conda_bcftools = os.path.join(sys.prefix, "bin", "bcftools")

    if os.path.isfile(conda_bcftools) and os.access(conda_bcftools, os.X_OK):
        tools["bcftools"] = conda_bcftools

    missing = [x for x, path in tools.items() if not path]
    if missing:
        raise FileNotFoundError("Missing independent branch tools: " + ", ".join(missing))
    threads = int(os.environ.get("AURORA_CMV_BCFTOOLS_THREADS", "4"))
    if not 1 <= threads <= 32:
        raise ValueError("AURORA_CMV_BCFTOOLS_THREADS must be 1..32")
    root = Path(sample_dir) / "cmv_bcftools_review"
    root.mkdir(parents=True, exist_ok=True)
    bam = root / (sample + ".merlin.sorted.bam")
    vcf = root / (sample + ".merlin.bcftools.norm.vcf.gz")
    tsv = root / (sample + ".merlin.bcftools.candidates.tsv")
    manifest = root / (sample + ".merlin.bcftools.manifest.json")
    log_path = root / (sample + ".bcftools.log")
    # Place a private copy of the reference in the review directory, rather
    # than writing an index into the user's reference library.
    local_ref = root / "NC_006273.2.merlin.fasta"
    shutil.copyfile(reference, local_ref)
    with log_path.open("w") as log:
        _checked([tools["samtools"], "faidx", local_ref], log)
        mapping = [tools["minimap2"], "-ax", "map-ont", "-t", threads, local_ref, fastq]
        sorting = [tools["samtools"], "sort", "-@", threads, "-o", bam, "-"]
        log.write("COMMAND: " + " ".join(map(str, mapping)) + " | " + " ".join(map(str, sorting)) + "\n")
        log.flush()
        mapper = subprocess.Popen(list(map(str, mapping)), stdout=subprocess.PIPE, stderr=log)
        try:
            subprocess.run(list(map(str, sorting)), stdin=mapper.stdout, check=True, stderr=log)
        finally:
            mapper.stdout.close()
        if mapper.wait() != 0:
            raise RuntimeError("minimap2 failed; inspect " + str(log_path))
        _checked([tools["samtools"], "index", bam], log)
        # Single-ploidy site calling: not a low-frequency minority-variant assay.
        pileup = [tools["bcftools"], "mpileup", "-Ou", "-f", local_ref,
                  "-q", "20", "-Q", "13", "-d", "100000", bam]
        calling = [tools["bcftools"], "call", "-m", "-v", "--ploidy", "1", "-Ou"]
        norm = [tools["bcftools"], "norm", "-f", local_ref, "-m", "-any", "-Oz", "-o", vcf]
        log.write("COMMAND: " + " | ".join(" ".join(map(str, x)) for x in (pileup, calling, norm)) + "\n")
        log.flush()
        p1 = subprocess.Popen(list(map(str, pileup)), stdout=subprocess.PIPE, stderr=log)
        p2 = subprocess.Popen(list(map(str, calling)), stdin=p1.stdout, stdout=subprocess.PIPE, stderr=log)
        p1.stdout.close()
        try:
            subprocess.run(list(map(str, norm)), stdin=p2.stdout, check=True, stderr=log)
        finally:
            p2.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
        if rc1 or rc2:
            raise RuntimeError(f"BCFtools pipeline failed: mpileup={rc1}, call={rc2}; inspect {log_path}")
        _checked([tools["bcftools"], "index", "-t", "-f", vcf], log)
        _checked([tools["bcftools"], "view", "-h", vcf], log, stdout=subprocess.DEVNULL)
        count = _vcf_to_tsv(tools["bcftools"], vcf, tsv, log)
    result = {
        "status": "REVIEW_ONLY", "sample": sample, "reference_accession": ACCESSION,
        "reference_source": str(reference), "original_fastq": str(fastq),
        "bam": str(bam), "normalized_vcf": str(vcf), "candidates_tsv": str(tsv),
        "log": str(log_path), "candidate_records": count,
        "limitations": [
            "Haploid consensus-oriented call, not validated for low-frequency mixtures or nanopore indels.",
            "Empty VCF is not evidence that resistance mutations are absent.",
            "Aurora protein calls are not directly comparable to VCF alleles without coordinate-level annotation.",
            "Independent remapping/calling only; shares the original FASTQ and is not orthogonal wet-lab validation.",
        ],
    }
    manifest.write_text(json.dumps(result, indent=2) + "\n")
    result["manifest"] = str(manifest)
    return result


def optional_run(sample, fastq, reference_dir, sample_dir):
    """Default-on for CMV, optional off via AURORA_CMV_BCFTOOLS=0."""
    if os.environ.get("AURORA_CMV_BCFTOOLS", "1").strip().lower() in ("0", "false", "no"):
        print(f"[CMV BCFTOOLS] {sample}: DISABLED (AURORA_CMV_BCFTOOLS=0)", flush=True)
        return {"cmv_bcftools_status": "DISABLED"}
    print(f"[CMV BCFTOOLS] {sample}: START independent Merlin mapping and BCFtools calling", flush=True)
    try:
        data = run(sample, fastq, reference_dir, sample_dir)
        print(f"[CMV BCFTOOLS] {sample}: REVIEW_ONLY; {data['candidate_records']} candidate records; VCF={data['normalized_vcf']}", flush=True)
        return {"cmv_bcftools_status": "REVIEW_ONLY", "cmv_bcftools_vcf": data["normalized_vcf"],
                "cmv_bcftools_manifest": data["manifest"], "cmv_bcftools_candidates": data["candidate_records"]}
    except Exception as exc:
        print(f"[CMV BCFTOOLS] {sample}: REVIEW_FAILED: {type(exc).__name__}: {exc}", flush=True)
        root = Path(sample_dir) / "cmv_bcftools_review"
        root.mkdir(parents=True, exist_ok=True)
        (root / (sample + ".bcftools.failure.json")).write_text(json.dumps({
            "status": "REVIEW_FAILED", "sample": sample, "error": str(exc),
            "error_type": type(exc).__name__}, indent=2) + "\n")
        return {"cmv_bcftools_status": "REVIEW_FAILED: " + str(exc)}
