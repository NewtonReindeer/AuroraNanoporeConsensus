# Aurora local BLAST setup on Ubuntu

Aurora now BLASTs each successfully generated Medaka consensus FASTA against a local NCBI nucleotide virus database during sample processing.

## Required command-line tools

The existing Aurora tools are still required, and `blastn` is now required as well:

```bash
which minimap2
which samtools
which medaka_consensus
which blastn
```

`blastdbcmd` is recommended because Aurora uses it to validate the database and record the database date/sequence count when available.

## Point Aurora at the database

Set the database **prefix**, not a `.nal`, `.nin`, or `.nsq` filename:

```bash
export AURORA_BLAST_DB=/path/to/AdenoPipeline/Database/nt_viruses
```

For example, if the Ubuntu system mounts NewtVault at `/mnt/NewtVault`:

```bash
export AURORA_BLAST_DB=/mnt/NewtVault/AdenoPipeline/Database/nt_viruses
```

Aurora also checks several common locations automatically, but the environment variable is the recommended production configuration.

Verify the database:

```bash
blastdbcmd -db "$AURORA_BLAST_DB" -info
```

## Performance controls

Optional environment variables:

```bash
export AURORA_BLAST_THREADS=8
export AURORA_BLAST_TIMEOUT_SECONDS=900
```

The default thread count is up to 8 CPUs and the default timeout is 900 seconds per sample.

## What is reported

For every virus, Aurora stores the BLAST top accession, inferred organism, identity, query coverage, E-value, bit score, database provenance, and the top-hit table.

For HBV and HCV, Aurora additionally derives a conservative genotype/subtype call from genotype-labelled high-quality top BLAST hits. This is reported separately from the mapping reference. HCV resistance-panel subtype selection remains based on the existing validated pipeline logic and is not silently changed by the BLAST call.

## Database updating

Database updating is intentionally not performed during a sequencing run. Keep the `nt_viruses` database maintenance process separate (for example with a systemd timer) so routine analyses do not block on a large download and each run has reproducible database provenance.
