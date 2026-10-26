"""FASTA sequence parsing and codon helpers."""
import os
import re


def read_single_fasta_sequence(fasta_path):
    seq_parts = []
    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith(">"):
                    seq_parts.append(line.strip())
        return "".join(seq_parts).upper()
    except Exception as e:
        print(f"Error reading FASTA sequence from {fasta_path}: {e}")
        return ""

def read_named_fasta_sequence(fasta_path, ref_name):
    """Read one named sequence from a possibly multi-FASTA reference file."""
    if not fasta_path or not os.path.exists(fasta_path) or not ref_name:
        return ""

    wanted = str(ref_name).split()[0]
    seq_parts = []
    collecting = False
    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    current = line[1:].strip().split()[0]
                    if collecting:
                        break
                    collecting = current == wanted
                    continue
                if collecting:
                    seq_parts.append(line.strip())
    except Exception as e:
        print(f"Error reading reference sequence {ref_name} from {fasta_path}: {e}")
        return ""

    return re.sub(r"[^ACGTN]", "", "".join(seq_parts).upper())

def fasta_sequence_stats(fasta_path):
    seq = read_single_fasta_sequence(fasta_path)
    if not seq:
        return {
            "consensus_length": 0,
            "n_count": 0,
            "n_percent": 100.0,
            "non_n_length": 0,
        }

    seq = re.sub(r"[^ACGTN]", "", seq.upper())
    total = len(seq)
    n_count = seq.count("N")
    non_n = total - n_count

    return {
        "consensus_length": total,
        "n_count": n_count,
        "n_percent": round((n_count / total) * 100, 2) if total else 100.0,
        "non_n_length": non_n,
    }

def amplicon_fasta_sequence_stats(fasta_path, start_1based, end_1based):
    """Calculate sequence/QC statistics across a detected amplicon interval.

    Unlike the generic FASTA helper, this function preserves standard IUPAC
    ambiguity characters. This matters when Aurora uses the BAM-derived IUPAC
    consensus as a safety fallback for a failed Medaka consensus.
    """
    seq = read_single_fasta_sequence(fasta_path)

    empty = {
        "consensus_length": 0,
        "n_count": 0,
        "n_percent": 100.0,
        "non_n_length": 0,
        "ambiguity_count": 0,
        "ambiguity_percent": 0.0,
        "acgt_count": 0,
        "acgt_percent": 0.0,
    }
    if not seq or not start_1based or not end_1based or end_1based < start_1based:
        return empty

    # BLAST and FASTA both support these standard nucleotide ambiguity codes.
    seq = re.sub(r"[^ACGTRYSWKMBDHVN]", "", seq.upper().replace("U", "T"))
    target_length = end_1based - start_1based + 1

    # If consensus retains full-reference coordinates, extract the amplicon.
    if len(seq) >= end_1based:
        target_seq = seq[start_1based - 1:end_1based]

    # If Medaka/IUPAC generation already emitted only the amplicon, use it directly.
    elif len(seq) <= target_length + 20:
        target_seq = seq

    else:
        target_seq = seq[:target_length]

    total = len(target_seq)
    n_count = target_seq.count("N")
    ambiguity_count = sum(target_seq.count(code) for code in "RYSWKMBDHV")
    acgt_count = sum(target_seq.count(base) for base in "ACGT")
    non_n = total - n_count

    return {
        "consensus_length": total,
        "n_count": n_count,
        "n_percent": round((n_count / total) * 100, 2) if total else 100.0,
        "non_n_length": non_n,
        "ambiguity_count": ambiguity_count,
        "ambiguity_percent": round((ambiguity_count / total) * 100, 2) if total else 0.0,
        "acgt_count": acgt_count,
        "acgt_percent": round((acgt_count / total) * 100, 2) if total else 0.0,
    }

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",

    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",

    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",

    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate_codon(codon):
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3:
        return "X"
    if any(base not in "ACGT" for base in codon):
        return "X"
    return CODON_TABLE.get(codon, "X")

def circular_subseq(seq, start_1based, length):
    if not seq:
        return ""
    n = len(seq)
    start_idx = (start_1based - 1) % n
    out = []
    for i in range(length):
        out.append(seq[(start_idx + i) % n])
    return "".join(out)

def circular_codon(seq, start_1based):
    return circular_subseq(seq, start_1based, 3)
