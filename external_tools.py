"""External executable discovery and optional subprocess integrations."""
import os
import subprocess
import sys

from .filesystem import app_base_dir, path_exists_any


def which_or_none(exe_name):
    try:
        import shutil
        return shutil.which(exe_name)
    except Exception:
        return None

def check_required_tools():
    missing = []
    for tool in ["minimap2", "samtools", "medaka_consensus", "blastn"]:
        if which_or_none(tool) is None:
            missing.append(tool)

    if missing:
        raise FileNotFoundError(
            "Missing required tools on PATH: " + ", ".join(missing) +
            "\nPlease install them and make sure they are available in your shell PATH."
        )

def safe_remove(path):
    try:
        if path and os.path.exists(path):
            os.remove(path)
    except Exception:
        pass

def run_ncbi_database_script(input_folder):
    try:
        base = app_base_dir()

        ncbi_script_path = path_exists_any(
            os.path.join(base, "AdenoPipeline", "BLASTscript.py"),
            "/Volumes/Newt-Vault/AdenoPipeline/BLASTscript.py",
            "/Volumes/NewtVault/AdenoPipeline/BLASTscript.py"
        )

        if not ncbi_script_path:
            print("NCBI script not found. Put BLASTscript.py in ./AdenoPipeline/ next to this script or use the mounted volume path.")
            return

        output_excel = os.path.join(input_folder, "BLAST_results.xlsx")
        invalid_files_csv = os.path.join(input_folder, "invalid_files.csv")

        blastplus_path = path_exists_any(
            os.path.join(base, "blastplus"),
            "/Users/xsjoak/Documents/Entero/blastplus"
        )

        if not blastplus_path:
            print("BLAST+ folder not found. Put it in ./blastplus next to the script or adjust the path.")
            return

        subprocess.run([
            sys.executable,
            ncbi_script_path,
            input_folder,
            blastplus_path,
            output_excel,
            invalid_files_csv
        ], check=True)

        print(f"NCBI database script executed successfully. Results saved to {output_excel}")

    except subprocess.CalledProcessError as e:
        print(f"Error running the NCBI database script: {e}")
    except Exception as e:
        print(f"Unexpected error while running the NCBI script: {e}")
