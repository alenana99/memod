# Gene Set Enrichment Analysis (GSEA)

This project performs gene set enrichment analysis using methylation data and GO annotations.

## Requirements

- R
- Python 3
- R packages: GO.db
- R packages: AnnotationDbi
- Python packages: pandas

## Instructions

1. Run `process_methylation.py` to count the number of methylations per gene.
2. Run `process_go.R` to process GO terms and create gene sets.
3. Run `run_gsea.R` to run the GSEA analysis and calculate ES, NES, and FDR.

## Execution

```bash
python scripts/process_methylation.py
Rscript scripts/process_go.R
Rscript scripts/run_gsea.R
