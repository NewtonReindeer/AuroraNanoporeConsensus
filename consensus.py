"""Medaka model selection and consensus generation."""
import gzip
import os
import re
import subprocess

from .cancellation import run_cancellable
from .filesystem import ensure_dir
from .fasta import read_named_fasta_sequence

try:
    import pysam
except ImportError:
    pysam = None


IUPAC_BASES = {
    frozenset({"A"}): "A",
    frozenset({"C"}): "C",
    frozenset({"G"}): "G",
    frozenset({"T"}): "T",
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",
    frozenset({"A", "C", "G", "T"}): "N",
}

def normalise_dorado_model_to_medaka(raw_model):
    """Convert common Dorado model identifiers to Medaka model names."""
    if not raw_model:
        return None
    model = str(raw_model).strip().strip('"').strip("'")
    model = model.split(',')[0].split(';')[0]
    model = re.sub(r"^(?:dna|rna)_", "", model, flags=re.IGNORECASE)
    model = model.replace("@", "_")
    model = re.sub(r"r10\.4\.1", "r1041", model, flags=re.IGNORECASE)
    model = re.sub(r"e8\.2", "e82", model, flags=re.IGNORECASE)
    model = re.sub(r"r9\.4\.1", "r941", model, flags=re.IGNORECASE)
    model = re.sub(r"[^A-Za-z0-9_.-]+", "_", model)
    return model.strip("_") or None

def detect_medaka_model_from_fastq(fastq_gz_path, max_records=50):
    """Inspect early FASTQ headers and return the most frequent detected model."""
    from collections import Counter
    detected = Counter()
    patterns = [
        re.compile(r"(?:model_version_id|basecall_model|basecaller_model|model)=([^\s]+)", re.I),
        re.compile(r"(?:model_version_id|basecall_model|basecaller_model|model):([^\s]+)", re.I),
    ]
    try:
        with gzip.open(fastq_gz_path, "rt", encoding="utf-8", errors="replace") as handle:
            for record_index in range(max_records):
                header = handle.readline()
                if not header:
                    break
                handle.readline(); handle.readline(); handle.readline()
                if not header.startswith("@"):
                    continue
                for pattern in patterns:
                    match = pattern.search(header)
                    if match:
                        normalised = normalise_dorado_model_to_medaka(match.group(1))
                        if normalised:
                            detected[normalised] += 1
                        break
    except Exception as exc:
        print(f"[MODEL WARNING] Could not inspect {fastq_gz_path}: {exc}")
        return None
    return detected.most_common(1)[0][0] if detected else None

def choose_medaka_model(input_fastq, fallback_model, auto_detect_models=False):
    if not auto_detect_models:
        return fallback_model
    detected = detect_medaka_model_from_fastq(input_fastq)
    if detected:
        print(f"[MODEL] Auto-detected Medaka model for {os.path.basename(input_fastq)}: {detected}")
        return detected
    print(
        f"[MODEL WARNING] No model metadata detected in {os.path.basename(input_fastq)}; "
        f"using fallback model: {fallback_model}"
    )
    return fallback_model

def generate_iub_consensus_from_bam(
    bam_file,
    ref_name,
    ref_path,
    output_folder,
    sample_name,
    variant_frequency=0.15,
    min_base_quality=20,
    min_mapping_quality=20,
    trim_start_1based=None,
    trim_end_1based=None,
):
    """Create a standalone IUPAC/IUB consensus FASTA from aligned reads.

    Bases present at or above variant_frequency are represented together with
    the standard IUPAC ambiguity code. Positions without usable read support
    are written as N. Insertions are not represented.

    For amplicon workflows, trim_start_1based/trim_end_1based restrict the
    exported IUB sequence to the detected read-supported target interval. This
    removes unsupported reference flanks while preserving Ns and IUPAC calls
    inside the actually sequenced region.
    """
    if pysam is None:
        print("[IUB WARNING] pysam is unavailable; IUB consensus was not created.")
        return None

    reference_sequence = read_named_fasta_sequence(ref_path, ref_name)
    if not reference_sequence:
        print(f"[IUB WARNING] Could not read reference '{ref_name}' from {ref_path}.")
        return None

    reference_sequence = reference_sequence.upper().replace("U", "T")
    calls = ["N"] * len(reference_sequence)

    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for column in bam.pileup(
                ref_name,
                0,
                len(reference_sequence),
                truncate=True,
                stepper="samtools",
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                max_depth=1000000,
            ):
                pos = column.reference_pos
                if pos < 0 or pos >= len(calls):
                    continue

                counts = {"A": 0, "C": 0, "G": 0, "T": 0}
                for pileup_read in column.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    alignment = pileup_read.alignment
                    if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                        continue
                    query_pos = pileup_read.query_position
                    if query_pos is None or not alignment.query_sequence:
                        continue
                    base = alignment.query_sequence[query_pos].upper().replace("U", "T")
                    if base in counts:
                        counts[base] += 1

                depth = sum(counts.values())
                if depth == 0:
                    continue

                represented = {
                    base for base, count in counts.items()
                    if count > 0 and (count / depth) >= variant_frequency
                }

                # Always retain the most abundant observed base if rounding or
                # an unusually mixed position leaves no base above threshold.
                if not represented:
                    represented = {max(counts, key=counts.get)}

                calls[pos] = IUPAC_BASES.get(frozenset(represented), "N")

        ensure_dir(output_folder, must_be_writable=True)
        output_path = os.path.join(output_folder, f"{sample_name}.consensus.iub.fasta")
        sequence = "".join(calls)

        trim_start = int(trim_start_1based) if trim_start_1based else None
        trim_end = int(trim_end_1based) if trim_end_1based else None
        trim_note = ""
        if trim_start and trim_end and trim_end >= trim_start:
            if trim_end <= len(sequence):
                sequence = sequence[trim_start - 1:trim_end]
                trim_note = f" amplicon={trim_start}-{trim_end}"
                print(
                    f"IUB consensus restricted to detected amplicon "
                    f"{trim_start}-{trim_end} ({len(sequence)} bases)."
                )
            else:
                print(
                    f"[IUB WARNING] Requested amplicon {trim_start}-{trim_end} "
                    f"does not fit reference-length IUB consensus ({len(sequence)} bases); "
                    "leaving sequence unchanged."
                )

        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write(
                f">{sample_name} IUB_consensus ref={ref_name} "
                f"threshold={variant_frequency:.3f}{trim_note}\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

        print(f"IUB consensus completed. Output: {output_path}")
        return output_path
    except Exception as exc:
        print(f"[IUB WARNING] Could not create IUB consensus for {sample_name}: {exc}")
        return None

def _read_fasta_records(fasta_path):
    """Read FASTA records as ``[(name, sequence), ...]`` without joining contigs."""
    records = []
    current_name = None
    seq_parts = []
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, "".join(seq_parts).upper()))
                current_name = line[1:].strip()
                seq_parts = []
            elif current_name is not None:
                seq_parts.append(line)
        if current_name is not None:
            records.append((current_name, "".join(seq_parts).upper()))
    return records


def generate_consensus(
    input_fastq,
    draft_fasta,
    medaka_model,
    output_folder,
    sample_name,
    ref_name=None,
    trim_start_1based=None,
    trim_end_1based=None,
):
    """Run Medaka and return one sample consensus FASTA.

    Medaka can emit more than one FASTA record when the draft/reference is a
    multi-FASTA. Aurora must never flatten those records into one artificial
    sequence. When ``ref_name`` is supplied, only the record matching the
    selected mapping reference is retained.

    For amplicon workflows, ``trim_start_1based``/``trim_end_1based`` restrict
    the exported sequence to the read-supported target interval. Internal Ns
    are deliberately preserved; only reference sequence outside the detected
    amplicon is excluded.
    """
    try:
        medaka_output_folder = os.path.join(output_folder, f"{sample_name}_medaka")
        ensure_dir(medaka_output_folder, must_be_writable=True)

        cmd = [
            "medaka_consensus",
            "-i", input_fastq,
            "-d", draft_fasta,
            "-o", medaka_output_folder,
            "-m", medaka_model,
        ]

        env = os.environ.copy()
        env["PYTHONNOUSERSITE"] = "1"
        run_cancellable(cmd, check=True, env=env)

        medaka_out = os.path.join(medaka_output_folder, "consensus.fasta")
        if not (os.path.exists(medaka_out) and os.path.getsize(medaka_out) > 0):
            print(f"Medaka did not produce the expected output file: {medaka_out}")
            return None

        records = _read_fasta_records(medaka_out)
        if not records:
            print(f"Warning: consensus for {sample_name} is empty.")
            return None

        selected_name = ""
        selected_seq = ""
        wanted = str(ref_name or "").split()[0]

        if wanted:
            for header, sequence in records:
                if header.split()[0] == wanted:
                    selected_name = header
                    selected_seq = sequence
                    break

        if not selected_seq:
            if len(records) == 1:
                selected_name, selected_seq = records[0]
            else:
                # Do not recreate the old behaviour of concatenating unrelated
                # Medaka contigs. The longest record is a safer fallback and is
                # clearly logged so unexpected reference naming is visible.
                selected_name, selected_seq = max(records, key=lambda item: len(item[1]))
                print(
                    f"[MEDAKA WARNING] Could not find selected reference '{wanted}' "
                    f"among {len(records)} Medaka contigs; using longest contig "
                    f"'{selected_name.split()[0]}' instead."
                )

        selected_seq = re.sub(r"[^ACGTN]", "", selected_seq.upper())
        if not selected_seq:
            print(f"Warning: selected Medaka consensus for {sample_name} contains no sequence.")
            return None

        trim_start = int(trim_start_1based) if trim_start_1based else None
        trim_end = int(trim_end_1based) if trim_end_1based else None
        if trim_start and trim_end and trim_end >= trim_start:
            target_length = trim_end - trim_start + 1
            if len(selected_seq) >= trim_end:
                selected_seq = selected_seq[trim_start - 1:trim_end]
                print(
                    f"Medaka consensus restricted to detected amplicon "
                    f"{trim_start}-{trim_end} ({len(selected_seq)} bases)."
                )
            elif len(selected_seq) <= target_length + 20:
                # Medaka has already emitted an amplicon-sized sequence.
                print(
                    f"Medaka consensus already appears amplicon-sized "
                    f"({len(selected_seq)} bases); no coordinate crop applied."
                )
            else:
                print(
                    f"[MEDAKA WARNING] Requested amplicon {trim_start}-{trim_end} "
                    f"does not fit selected consensus length {len(selected_seq)}; "
                    "leaving sequence unchanged."
                )

        renamed = os.path.join(medaka_output_folder, f"{sample_name}.consensus.fasta")
        with open(renamed, "w", encoding="utf-8") as handle:
            ref_note = f" ref={wanted}" if wanted else ""
            handle.write(f">{sample_name}{ref_note}\n")
            for seq_start in range(0, len(selected_seq), 80):
                handle.write(selected_seq[seq_start:seq_start + 80] + "\n")

        print(f"Medaka completed. Output: {renamed}")
        return renamed if os.path.getsize(renamed) > 0 else None

    except subprocess.CalledProcessError as e:
        print(f"Error in Medaka process for {input_fastq}: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error in Medaka for {input_fastq}: {e}")
        return None
