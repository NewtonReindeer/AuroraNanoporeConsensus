"""Static and smoke validation for the first-stage Aurora module split."""
from __future__ import annotations

import ast
import compileall
import importlib
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
ORIGINAL = PACKAGE_DIR / "original_monolith.py"
EXCLUDED = {"original_monolith.py", "legacy_api.py", "app.py", "__main__.py", "validate_refactor.py"}

# Definitions intentionally changed by the integrated local-BLAST feature.
# They are still required to exist exactly once, but are no longer expected to
# be AST-identical to the audit monolith.
INTENTIONALLY_CHANGED = {
    "check_required_tools",
    "process_sample",
    "run_full_pipeline",
    "log_reference_usage",
    "log_hcv_resistance_results",
    "log_hbv_mutation_results",
    "write_hcv_html_report",
    "generate_hcv_html_reports",
    "write_hbv_html_report",
    "AuroraWindow",
}


def top_level_definitions(path: Path) -> dict[str, ast.AST]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
    }


def validate_definition_mapping() -> None:
    original = top_level_definitions(ORIGINAL)
    original.pop("main", None)  # app.py intentionally uses a lazy-GUI equivalent.

    moved: dict[str, tuple[Path, ast.AST]] = {}
    duplicates: dict[str, list[str]] = {}
    for path in PACKAGE_DIR.glob("*.py"):
        if path.name in EXCLUDED:
            continue
        for name, node in top_level_definitions(path).items():
            if name in original and name in moved:
                duplicates.setdefault(name, [moved[name][0].name]).append(path.name)
            moved[name] = (path, node)

    missing = sorted(set(original) - set(moved))
    if missing or duplicates:
        raise AssertionError(f"Missing={missing}; duplicates={duplicates}")

    changed = []
    for name, original_node in original.items():
        if name in INTENTIONALLY_CHANGED:
            continue
        moved_node = moved[name][1]
        if ast.dump(original_node, include_attributes=False) != ast.dump(moved_node, include_attributes=False):
            changed.append(name)
    if changed:
        raise AssertionError(f"Moved definitions changed structurally: {changed}")

    print(f"Definition mapping: unchanged moved definitions are AST-identical; intentional changes: {sorted(INTENTIONALLY_CHANGED)}")


def validate_compilation() -> None:
    if not compileall.compile_dir(str(PACKAGE_DIR), quiet=1):
        raise AssertionError("Compilation failed")
    print("Compilation: OK")


def validate_non_gui_imports() -> None:
    modules = [
        "alignment", "blast", "config", "consensus", "external_tools", "fasta",
        "filesystem", "hbv", "hcv", "qc", "reporting_html",
        "reporting_excel", "pipeline", "legacy_api",
    ]
    for module in modules:
        importlib.import_module(f"{__package__}.{module}")
    print("Non-GUI imports: OK")


def main() -> int:
    validate_compilation()
    validate_definition_mapping()
    validate_non_gui_imports()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
