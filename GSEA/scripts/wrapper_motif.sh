#!/bin/bash

INSTALLATION_DIR=/home/alessiam/memod/Source
circ="$INSTALLATION_DIR/circular_plotter.R"

# Controlla il numero di argomenti passati
if [ $# -lt 3 ]; then
    echo "Usage: $0 <motif_dir1> <motif_dir2> ... <motif_dirN> <output_dir>"
    exit 1
fi

# Percorsi delle directory dei motivi (i primi argomenti)
motif_dirs=("${@:1:$#-1}")

# Directory di output (l'ultimo argomento)
output_dir="${@: -1}"

# Ciclo su ogni directory dei motivi
for motif_dir in "${motif_dirs[@]}"
do
    # Estrai il nome del motivo dalla directory
    motif_name=$(basename "$motif_dir")

    # File BED per il motivo corrente
    bed_files=(
        "${motif_dir}/msa/${motif_name}_CDS.bed"
        "${motif_dir}/msa/${motif_name}_nCDS.bed"
        "${motif_dir}/msa/${motif_name}_tIG.bed"
        "${motif_dir}/msa/${motif_name}_US.bed"
    )

  # Chiamata allo script R usando Rscript
      Rscript "$circ" "${bed_files[0]}" "${bed_files[1]}" "${bed_files[2]}" "${bed_files[3]}"

done
