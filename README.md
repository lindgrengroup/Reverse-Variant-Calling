# Reverse Variant Calling
Workflow for SNP calling from GTEX RNA-seq data using the GATK Best Practices pipeline.

The RNA-seq data obtained (.bam files, not provided here) were from the GTEX project, V8. They had been aligned to the reference genome GRCh38 using STAR v2.5.3a.

The industry-standard GATK Best Practices was closely followed (with the help of the UCLA workshop) with the addition of SplitNCigar for splitting alignment overlapping exon/intron junctions and rescaling mapping quality.

## Resources:
UCLA youtube tutorial 1: https://www.youtube.com/watch?v=xLm3Le0rfYQ <br />
UCLA youtube tutorial 2: https://www.youtube.com/watch?v=hRsjy1Z8QDA&t=552s <br />
GATK pipeline: https://gatk.broadinstitute.org/hc/en-us <br />
RNA-calling: https://link.springer.com/protocol/10.1007/978-1-0716-2293-3_13 <br />

## Other software packages used:
- https://github.com/samtools/
- https://broadinstitute.github.io/picard/
