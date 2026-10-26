# CMV pipeline — initial implementation

The CMV selector routes through Aurora's existing reference selection, alignment, consensus, BLAST, and QC pipeline. A CMV-specific module detects exploratory minor SNVs (>=15% allele fraction, >=30x depth, minimum base quality 20 and mapping quality 20). Reports are HTML and XLSX.

**Not yet clinically validated:** Gene mapping, amino-acid numbering, reference-coordinate verified resistance calls, indels, haplotypes, strand-bias filtering, replicate support, and coverage-per-resistance-codon remain to be implemented. An optional `cmv_mutation_panel.tsv` can be parsed, but is **not used for clinical annotation** until coordinate mapping and validation are implemented. No drug susceptibility inference is made.

The detected target region and reference breadth are reported independently. The pipeline must not claim whole-genome recovery from an amplicon. CMV BLAST database setup, FASTA selection, external dependencies, and actual sequencing data require integration testing on the deployment machine.

Relevant genes to curate: UL97, UL54, UL56, UL89 and UL51. Use curated, versioned clinical evidence and exact reference accessions; do not treat an arbitrary genomic SNV as a drug-resistance mutation.
