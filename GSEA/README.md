# Gene Set Enrichment Analysis (GSEA)

This project performs gene set enrichment analysis using methylation data and GO annotations.

## Requirements

- R
- Python 3
- R packages: GO.db
- Python packages: pandas

## Instructions

1. Run `fetch_go_details.R` to get details about GO terms.
2. Run `count_methylations.py` to count the number of methylations per gene.
3. Run `process_go_terms.R` to process GO terms and create gene sets.
4. Run `create_gene_sets.R` to create the final gene sets for GSEA.
5. Run `run_gsea.R` to run the GSEA analysis and calculate ES, NES, and FDR.

## Execution

```bash
Rscript scripts/fetch_go_details.R
python scripts/count_methylations.py
Rscript scripts/process_go_terms.R
Rscript scripts/create_gene_sets.R
Rscript scripts/run_gsea.R
