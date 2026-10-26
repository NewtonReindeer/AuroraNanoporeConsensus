# Aurora refactor

This directory contains a behaviour-preserving first-stage split of the original
4,934-line Aurora script. The original source is retained as
`original_monolith.py` for audit and rollback.

## Run

Place the generated Qt files and resources beside these modules:

- `ui_aurora_mainwindow.py`
- `resources_rc.py`
- any referenced images/metadata files

Then run from the directory containing the `aurora_refactored` package:

```bash
python -m aurora_refactored
```

Required external programs remain unchanged: `minimap2`, `samtools`, and
`medaka_consensus`. Python dependencies remain unchanged (`PySide6`, `openpyxl`,
and optional/required-for-variant-calling `pysam`).

## Module map

- `app.py`, `gui.py`: application startup and Qt code
- `config.py`: unchanged thresholds, panels, choices, and report constants
- `models.py`: documentation-only structural types; runtime dictionaries remain unchanged
- `filesystem.py`: paths, discovery, run identity, exports, FASTA file operations
- `external_tools.py`: executable checks and optional NCBI subprocess
- `alignment.py`: SAM/BAM parsing and depth calculations
- `fasta.py`: sequence and codon helpers
- `consensus.py`: Medaka and IUPAC consensus generation
- `qc.py`: general and HBV system QC
- `hcv.py`, `hbv.py`: virus-specific scientific analysis
- `reporting_excel.py`, `reporting_html.py`: output generation
- `pipeline.py`: sample processing and run orchestration
- `legacy_api.py`: broad re-export facade for migration

## Deliberately unchanged

Clinical thresholds, mutation panels, output field names, subprocess arguments,
report wording, GUI workflow, and per-sample pipeline order were copied without
scientific redesign. The only structural changes are imports and module location.
Because functions now live in separate modules, diagnostic strings that print
`__file__` identify the responsible module rather than the former monolithic file.

## Validate the split

From the directory containing the package:

```bash
python -m aurora_refactored.validate_refactor
```

This compiles every module, confirms that all 113 moved top-level definitions are
present exactly once and AST-identical to the audit copy, and imports every
non-GUI module.
