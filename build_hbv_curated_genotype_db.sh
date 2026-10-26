#!/usr/bin/env bash
set -euo pipefail

# Build the Aurora curated HBV genotype BLAST database from the existing local
# HBV small database. This does not download anything from the internet.
#
# Usage:
#   ./build_hbv_curated_genotype_db.sh /path/to/aurora_refactored
#
# Optional environment override:
#   AURORA_HBV_SOURCE_DB=/path/to/hbv/database/prefix

if [ "$#" -gt 1 ]; then
  echo "Usage: $0 [/path/to/aurora_refactored]" >&2
  exit 2
fi

if [ "$#" -eq 1 ]; then
  AURORA_DIR="$(cd "$1" && pwd)"
else
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  if [ -d "$SCRIPT_DIR/blastdb" ]; then
    AURORA_DIR="$SCRIPT_DIR"
  elif [ -d "$SCRIPT_DIR/../aurora_refactored" ]; then
    AURORA_DIR="$(cd "$SCRIPT_DIR/../aurora_refactored" && pwd)"
  else
    echo "Could not auto-detect aurora_refactored. Pass it as the first argument." >&2
    exit 2
  fi
fi

BLASTDBCMD="$(command -v blastdbcmd || true)"
MAKEBLASTDB="$(command -v makeblastdb || true)"
if [ -z "$BLASTDBCMD" ] || [ -z "$MAKEBLASTDB" ]; then
  echo "blastdbcmd and makeblastdb must be available on PATH." >&2
  exit 1
fi

SOURCE_DB="${AURORA_HBV_SOURCE_DB:-$AURORA_DIR/blastdb/aurora_small_db/hbv/hbv}"
OUT_DIR="$AURORA_DIR/blastdb/hbv_genotype_db"
OUT_PREFIX="$OUT_DIR/hbv_genotypes"
FASTA="$OUT_DIR/hbv_genotypes.fasta"
MANIFEST="$OUT_DIR/hbv_genotype_manifest.tsv"
TMP_FASTA="$OUT_DIR/.hbv_genotypes.tmp.fasta"

if ! "$BLASTDBCMD" -db "$SOURCE_DB" -info >/dev/null 2>&1; then
  echo "Cannot open source HBV BLAST database: $SOURCE_DB" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
: > "$TMP_FASTA"
printf 'genotype\tsubgenotype\taccession\n' > "$MANIFEST"

# Published full-genome human HBV genotype/subgenotype reference set.
# The classifier compares the best representative per genotype, so genotypes
# with more subgenotypes do not gain a voting advantage.
REFS=(
  'A|A1|JN182318'
  'A|A2|HE576989'
  'A|A3|AB194951'
  'A|A4|AY934764'
  'A|A5|FJ692613'
  'A|A6|GQ331047'
  'A|A7|FN545833'
  'B|B1|AB642091'
  'B|B2|FJ899779'
  'B|B3|GQ924617'
  'B|B4|GQ924626'
  'B|B5|GQ924640'
  'B|B6|JN792893'
  'B|B7|GQ358137'
  'B|B8|GQ358147'
  'B|B9|GQ358149'
  'C|C1|AB697490'
  'C|C2|GQ358158'
  'C|C3|DQ089801'
  'C|C4|HM011493'
  'C|C5|EU410080'
  'C|C6|EU670263'
  'C|C7|GU721029'
  'C|C8|AP011106'
  'C|C9|AP011108'
  'C|C10|AB540583'
  'C|C11|AB554019'
  'C|C12|AB554025'
  'C|C13|AB644280'
  'C|C14|AB644284'
  'C|C15|AB644286'
  'C|C16|AB644287'
  'D|D1|GU456636'
  'D|D2|GQ477452'
  'D|D3|EU594434'
  'D|D4|GQ922003'
  'D|D5|GQ205377'
  'D|D6|KF170740'
  'D|D7|FJ904442'
  'D|D8|FN594770'
  'D|D9|JN664942'
  'E|E|FN594748'
  'F|F1B|FJ709464'
  'F|F2B|DQ899146'
  'F|F1A|AY090459'
  'F|F2A|DQ899142'
  'F|F3|AB036920'
  'F|F4|AF223965'
  'G|G|GU563556'
  'H|H|AB516393'
  'I|I1|FJ023659'
  'I|I2|FJ023664'
  'J|J|AB486012'
)

missing=0
count=0
for spec in "${REFS[@]}"; do
  IFS='|' read -r genotype subgenotype accession <<< "$spec"
  echo "Extracting $accession (genotype $genotype, $subgenotype)"

  sequence="$($BLASTDBCMD -db "$SOURCE_DB" -entry "$accession" -outfmt '%s' 2>/dev/null | tr -d '[:space:]' || true)"
  if [ -z "$sequence" ]; then
    # Some databases preserve a versioned accession only. Try accession.1 as a
    # pragmatic fallback; if that fails, abort rather than silently building an
    # incomplete clinical-style reference set.
    sequence="$($BLASTDBCMD -db "$SOURCE_DB" -entry "${accession}.1" -outfmt '%s' 2>/dev/null | tr -d '[:space:]' || true)"
  fi

  if [ -z "$sequence" ]; then
    echo "  MISSING: $accession" >&2
    missing=$((missing + 1))
    continue
  fi

  safe_sub="$(printf '%s' "$subgenotype" | tr -cd 'A-Za-z0-9_.-')"
  safe_acc="$(printf '%s' "$accession" | tr -cd 'A-Za-z0-9_.-')"
  printf '>AURORA_HBV_%s_%s_%s genotype=%s subgenotype=%s accession=%s\n' \
    "$genotype" "$safe_sub" "$safe_acc" "$genotype" "$subgenotype" "$accession" >> "$TMP_FASTA"
  printf '%s\n' "$sequence" | fold -w 80 >> "$TMP_FASTA"
  printf '%s\t%s\t%s\n' "$genotype" "$subgenotype" "$accession" >> "$MANIFEST"
  count=$((count + 1))
done

if [ "$missing" -ne 0 ]; then
  echo >&2
  echo "ERROR: $missing curated HBV reference accession(s) were not found in:" >&2
  echo "  $SOURCE_DB" >&2
  echo "No curated database was built. The partial FASTA is left at:" >&2
  echo "  $TMP_FASTA" >&2
  exit 1
fi

mv "$TMP_FASTA" "$FASTA"
rm -f "$OUT_PREFIX".n* "$OUT_PREFIX".p* 2>/dev/null || true

"$MAKEBLASTDB" \
  -in "$FASTA" \
  -dbtype nucl \
  -blastdb_version 5 \
  -title 'Aurora curated HBV genotype references v1' \
  -out "$OUT_PREFIX"

echo
echo "Curated HBV genotype database built successfully."
echo "Aurora directory: $AURORA_DIR"
echo "Source database:  $SOURCE_DB"
echo "Output prefix:    $OUT_PREFIX"
echo "References:       $count"
"$BLASTDBCMD" -db "$OUT_PREFIX" -info
