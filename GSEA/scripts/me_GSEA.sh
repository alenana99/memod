#!/bin/bash

# Directory contenente gli script e i binari
INSTALLATION_DIR="/home/alessiam/memod/prova_script"
process_go="$INSTALLATION_DIR/process_GO.R"
process_meth="$INSTALLATION_DIR/process_meth.py"
ES_wrapper_fast="$INSTALLATION_DIR/ES_wrapper_fast"
ES_significance2="$INSTALLATION_DIR/ES_significance2"
GSEA_frd="$INSTALLATION_DIR/GSEA_fdr"
FDR_wrapper="$INSTALLATION_DIR/FDR_wapper"
run_GSEA="$INSTALLATION_DIR/run_GSEA.R"

# Funzione per stampare l'help
print_help() {
  echo "Usage: me_gsea -go <str> -anno <str> -m <str> -o <str>"
  echo
  echo "-go    <str>           output PANNZER file"
  echo "-anno  <str>           prokka annotation file"
  echo "-m     <str>           tsv methylation file from MicrobeMod"
  echo "-o     <str>           output directory"
  echo
}

# Variabili per gli argomenti
GO_OUT=""
GENOMIC_GFF=""
METHYLATION_TSV=""
OUTPUT_DIR=""

# Analisi degli argomenti della riga di comando
while [[ $# -gt 0 ]]; do
  case $1 in
    -go)
      GO_OUT=$(readlink -m "$2")
      shift
      shift
      ;;
    -anno)
      GENOMIC_GFF=$(readlink -m "$2")
      shift
      shift
      ;;
    -m)
      METHYLATION_TSV=$(readlink -m "$2")
      shift
      shift
      ;;
    -o)
      OUTPUT_DIR=$(readlink -m "$2")
      shift
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      print_help
      exit 1
      ;;
    *)
      shift # past argument
      ;;
  esac
done

# Verifica degli argomenti richiesti
if [[ -z "$GO_OUT" || -z "$GENOMIC_GFF" || -z "$METHYLATION_TSV" || -z "$OUTPUT_DIR" ]]; then
  print_help
  exit 1
fi

# Esecuzione dello script Python per processare i dati di metilazione
python3 "$process_meth" "$GENOMIC_GFF" "$METHYLATION_TSV" "$OUTPUT_DIR/L"

# Esecuzione dello script R per processare i dati GO
Rscript "$process_go" "$GO_OUT" "$OUTPUT_DIR/methylation_output.csv" "$OUTPUT_DIR/list_S"

# Impostazione delle variabili base per i nomi dei file
base_go_out=$(basename "$GO_OUT")
base_methylation_tsv=$(basename "$METHYLATION_TSV")
base_genomic_gff=$(basename "$GENOMIC_GFF")

# Cambio della directory di output
cd "$OUTPUT_DIR" || exit

# Aggiornamento dei percorsi relativi nella directory di output
GENOMIC_GFF=$(readlink -m "$base_genomic_gff")
METHYLATION_TSV=$(readlink -m "$base_methylation_tsv")
GO_OUT=$(readlink -m "$base_go_out")

echo "Process completed. Output files are located in $OUTPUT_DIR"
