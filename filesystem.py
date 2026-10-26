"""Filesystem, run identity, and file-export helpers."""
import base64
import gzip
import mimetypes
import os
import re
import shutil
from datetime import datetime

from .config import REPORT_LOGO_FILENAMES


def app_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def ensure_dir(dirpath, must_be_writable=False):
    os.makedirs(dirpath, exist_ok=True)
    if must_be_writable:
        try:
            testfile = os.path.join(dirpath, ".write_test.tmp")
            with open(testfile, "w", encoding="utf-8") as f:
                f.write("ok")
            os.remove(testfile)
        except Exception as e:
            raise PermissionError(f"Cannot write to '{dirpath}': {e}")

def path_exists_any(*paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None

def is_real_fastq_gz(path):
    name = os.path.basename(path)
    if name.startswith("._") or name.startswith("."):
        return False
    return name.endswith(".fastq.gz")

def is_real_fasta(path):
    name = os.path.basename(path)
    if name.startswith("._") or name.startswith("."):
        return False
    return name.endswith(".fa") or name.endswith(".fasta")

def find_files_recursively(directory, predicate):
    matched = []
    try:
        for root, _, files in os.walk(directory):
            for file in files:
                full = os.path.join(root, file)
                try:
                    if predicate(full):
                        matched.append(full)
                except Exception:
                    continue
    except Exception as e:
        print(f"Error finding files in {directory}: {e}")
    return matched

def count_fastq_bases(fastq_gz_path):
    total_bases = 0
    try:
        with gzip.open(fastq_gz_path, "rt", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total_bases += len(line.strip())
    except Exception as e:
        print(f"Error counting FASTQ bases in {fastq_gz_path}: {e}")
    return total_bases

def build_run_tag(run_date, run_number):
    return f"{run_date}_run{int(run_number):02d}"

def allocate_run_identity(output_folder):
    """Allocate one date/run identifier shared by every output from this launch."""
    ensure_dir(output_folder, must_be_writable=True)
    run_date = datetime.now().strftime("%Y-%m-%d")
    pattern = re.compile(rf"^Aurora_{re.escape(run_date)}_run(\d+)\.run$")
    highest = 0
    for name in os.listdir(output_folder):
        match = pattern.match(name)
        if match:
            highest = max(highest, int(match.group(1)))
    run_number = highest + 1
    run_tag = build_run_tag(run_date, run_number)
    marker_path = os.path.join(output_folder, f"Aurora_{run_tag}.run")
    with open(marker_path, "x", encoding="utf-8") as marker_file:
        marker_file.write(
            f"date={run_date}\nrun_number={run_number}\n"
            f"created={datetime.now().isoformat(timespec='seconds')}\n"
            f"script={os.path.abspath(__file__)}\n"
        )
    return run_date, run_number, run_tag, marker_path

def safe_report_filename(sample):
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(sample or "prov")).strip("._")
    return safe or "prov"

def run_samples_folder(output_folder, run_tag):
    """Return the shared container for all per-sample/run subfolders."""
    return os.path.join(output_folder, f"Samples_{run_tag}")

def run_consensus_folder(output_folder, run_tag):
    """Return the run-level folder containing easy-to-find consensus copies."""
    return os.path.join(output_folder, f"Consensus_{run_tag}")

def is_valid_report_image(path):
    """Accept real PNG/JPEG files and reject macOS AppleDouble/resource-fork files."""
    if not path or not os.path.isfile(path):
        return False
    name = os.path.basename(path)
    if name.startswith("._") or name.startswith("."):
        return False
    try:
        with open(path, "rb") as handle:
            signature = handle.read(12)
        return (
            signature.startswith(b"\x89PNG\r\n\x1a\n")
            or signature.startswith(b"\xff\xd8\xff")
        )
    except Exception:
        return False

def export_run_consensus_files(result_rows, output_folder, run_tag):
    """Copy each sample's final consensus files into one run-level folder."""
    consensus_dir = run_consensus_folder(output_folder, run_tag)
    ensure_dir(consensus_dir, must_be_writable=True)
    copied = 0
    for row in result_rows:
        sample = safe_report_filename(row.get("sample", "sample"))
        sources = [
            (row.get("consensus_file"), f"{sample}.consensus.fasta"),
            (row.get("iub_consensus_file"), f"{sample}.consensus.iub.fasta"),
        ]
        for source, filename in sources:
            if not source or not os.path.isfile(source) or os.path.getsize(source) <= 0:
                continue
            destination = os.path.join(consensus_dir, filename)
            try:
                shutil.copy2(source, destination)
                copied += 1
            except Exception as exc:
                print(f"Could not copy consensus file {source} to {destination}: {exc}")
    if copied == 0:
        try:
            os.rmdir(consensus_dir)
        except OSError:
            pass
        print("No per-sample consensus files were available for run-level export.")
        return ""
    print(f"Per-sample consensus files saved to {consensus_dir}")
    return consensus_dir

def find_hbv_report_logo(row, report_path):
    """Find a PNG/JPEG logo near the reference, script, or report output."""
    candidates = []
    explicit = str(row.get("report_logo", "") or "").strip()
    if explicit:
        candidates.append(explicit)
    directories = [
        os.path.dirname(str(row.get("reference_file", "") or "")),
        app_base_dir(),
        os.path.dirname(report_path),
        os.path.dirname(os.path.dirname(report_path)),
    ]
    for directory in directories:
        if not directory:
            continue
        for filename in REPORT_LOGO_FILENAMES:
            candidates.append(os.path.join(directory, filename))
    for path in candidates:
        if is_valid_report_image(path):
            return path
    # Fallback: accept a clearly named logo PNG, or the only PNG in the reference folder.
    for directory in directories:
        if not directory or not os.path.isdir(directory):
            continue
        pngs = sorted(
            os.path.join(directory, name)
            for name in os.listdir(directory)
            if name.lower().endswith((".png", ".jpg", ".jpeg"))
            and not name.startswith("._")
            and not name.startswith(".")
            and is_valid_report_image(os.path.join(directory, name))
        )
        named_logos = [path for path in pngs if "logo" in os.path.basename(path).lower()]
        if named_logos:
            return named_logos[0]
        if len(pngs) == 1:
            return pngs[0]
    return ""

def image_file_to_data_uri(path):
    if not path:
        return ""
    try:
        mime = mimetypes.guess_type(path)[0] or "image/png"
        with open(path, "rb") as handle:
            encoded = base64.b64encode(handle.read()).decode("ascii")
        return f"data:{mime};base64,{encoded}"
    except Exception as exc:
        print(f"Kunde inte läsa rapportlogotyp {path}: {exc}")
        return ""

def concatenate_fasta_files(fasta_files, output_file):
    if not fasta_files:
        print("No consensus FASTA files to concatenate.")
        return False

    try:
        ensure_dir(os.path.dirname(output_file), must_be_writable=True)
        with open(output_file, "w", encoding="utf-8") as outfile:
            for fasta in fasta_files:
                if os.path.isfile(fasta) and os.path.getsize(fasta) > 0:
                    with open(fasta, "r", encoding="utf-8", errors="replace") as infile:
                        outfile.write(infile.read())
                else:
                    print(f"Skipping empty or missing FASTA file: {fasta}")

        print(f"Concatenated FASTA file saved to {output_file}")
        return True

    except Exception as e:
        print(f"Error concatenating FASTA files: {e}")
        return False

def split_fasta(output_fasta, split_folder):
    ensure_dir(split_folder, must_be_writable=True)

    with open(output_fasta, "r", encoding="utf-8", errors="replace") as f:
        current = None
        for line in f:
            if line.startswith(">"):
                if current:
                    current.close()
                header = line.strip()[1:].replace(" ", "_")
                current = open(os.path.join(split_folder, f"{header}.fasta"), "w", encoding="utf-8")
            if current:
                current.write(line)

        if current:
            current.close()

    print(f"Split FASTA sequences saved to {split_folder}")
