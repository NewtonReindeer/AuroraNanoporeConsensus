# Aurora — Nanopore Viral Sequencing Analysis

**Aurora** is a Python-based bioinformatics application designed to support viral genomic analysis using Oxford Nanopore sequencing data.

Aurora provides a graphical user interface for processing sequencing data, generating consensus sequences, performing quality control, and producing analysis reports.

The software includes dedicated analysis workflows for hepatitis B virus (HBV), hepatitis C virus (HCV), and cytomegalovirus (CMV).

## Features

- **Graphical user interface:** Desktop application built using PySide6 (Qt).
- **Sequencing analysis:** Processing of viral sequencing data using established bioinformatics tools.
- **Consensus generation:** Reference-based sequence alignment and consensus generation using Minimap2, SAMtools, and Medaka.
- **Quality control:** Assessment of sequencing quality, coverage, and analysis metrics.
- **HBV analysis:** HBV-specific analysis and genotyping functionality.
- **HCV analysis:** HCV-specific sequence analysis.
- **CMV analysis:** Dedicated CMV analysis modules, including variant calling, sequence comparison, and mutation annotation.
- **Reporting:** Generation of HTML and Excel analysis reports.
- **Local reference databases:** Support for NCBI BLAST databases and virus-specific reference datasets.

## System requirements

Aurora is a Python application that relies on several external bioinformatics programs.

### Python dependencies

- Python 3 (specific supported versions should be confirmed)
- PySide6
- openpyxl
- pysam (required for relevant variant-calling workflows)

Install Python dependencies using:

```bash
pip install -r requirements.txt
```

### External bioinformatics tools

The following programs are required for the relevant analysis workflows:

- Minimap2 — sequence alignment
- SAMtools — alignment file processing
- Medaka (`medaka_consensus`) — consensus generation
- NCBI BLAST+ — database searching and reference-based analysis where used

These tools must be installed separately and accessible through the system PATH, or configured according to Aurora's executable settings.

Additional workflow-specific dependencies may be required. See `BLAST_SETUP_UBUNTU.md` for BLAST-related setup information.

## Installation

Clone the repository:

```bash
git clone https://github.com/NewtonReindeer/AuroraNanoporeConsensus.git
```

Enter the repository directory:

```bash
cd AuroraNanoporeConsensus
```

Create a Python virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

### Launching Aurora

Aurora is organized as a Python package.

From the parent directory of the `aurora_refactored` package, run:

```bash
python -m aurora_refactored
```

The package must include the generated Qt interface and resource files, including `ui_aurora_mainwindow.py` and `resources_rc.py`.

## BLAST reference databases

**Large reference databases are not included in this repository.**

Aurora uses local BLAST databases and reference datasets for virus-specific analysis. The complete database installation can require approximately 70–80 GB of storage.

To keep the GitHub repository manageable, the large database files have been excluded.

The expected local directory structure includes:

```text
blastdb/
├── aurora_small_db/
├── hbv_genotype_db/
├── cmv_genotype_db/
├── cmv_genotype_db_v3/
└── nt_viruses.*
```

Users must obtain and prepare the relevant reference databases separately.

NCBI provides BLAST database downloads and documentation:

https://ftp.ncbi.nlm.nih.gov/blast/db/

https://www.ncbi.nlm.nih.gov/books/NBK279690/

Some Aurora-specific curated reference databases may require additional preparation using the scripts supplied with the project.

See `BLAST_SETUP_UBUNTU.md` and `build_hbv_curated_genotype_db.sh` for additional information.

**Important:** Downloading NCBI databases alone may not reproduce all of Aurora's curated reference datasets. Exact dataset versions and preparation steps should be documented for reproducibility.

## Project structure

The application is organized into modules covering different parts of the sequencing workflow.

| Module | Function |
|---|---|
| `app.py`, `gui.py` | Application startup and graphical interface |
| `pipeline.py` | Analysis workflow orchestration |
| `alignment.py` | Sequence alignment processing and coverage analysis |
| `consensus.py` | Consensus sequence generation |
| `qc.py` | Quality control |
| `hbv.py`, `hbv_genotyping.py` | HBV analysis and genotyping |
| `hcv.py` | HCV-specific analysis |
| `cmv.py`, `cmv_workflow.py` | CMV analysis workflows |
| `reporting_html.py` | HTML report generation |
| `reporting_excel.py` | Excel report generation |
| `config.py` | Analysis parameters and configuration |

Additional modules provide virus-specific analysis, mutation annotation, validation, and supporting functionality.

## Validation

The repository includes scripts for checking the refactored codebase.

Run the refactor validation utility from the package's parent directory:

```bash
python -m aurora_refactored.validate_refactor
```

Additional CMV-related test scripts are included in the repository.

These checks do not, by themselves, establish clinical validation of the complete application.

## Intended use

Aurora is developed for viral sequencing analysis in a laboratory bioinformatics environment.

**Clinical use requires appropriate local validation, quality assurance, and review of results.** Outputs should be interpreted in the context of validated laboratory workflows.

## Development status

Aurora is under active development. This repository contains the refactored Aurora 3.2 codebase, including HBV, HCV, and CMV analysis functionality.

The original monolithic implementation is retained in `original_monolith.py` for audit and comparison.

## Repository

Developed and maintained by [NewtonReindeer](https://github.com/NewtonReindeer).

Repository: https://github.com/NewtonReindeer/AuroraNanoporeConsensus
