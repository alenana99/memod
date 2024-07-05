#!/bin/bash

# INSTALLATION_DIR holds the path to the scripts, make sure you add the installation directory to your PATH

INSTALLATION_DIR=/home/alessiam/memod/src
#INSTALLATION_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
circ="$INSTALLATION_DIR/circular_plotter.R"

# Check the number of arguments passed
if [ $# -lt 3 ]; then
    echo "Usage: $0 <motif_dir1> <motif_dir2> ... <motif_dirN> <output_dir>"
    exit 1
fi

# Motif directory paths (the first arguments)
motif_dirs=("${@:1:$#-1}")

# Output Directory (the last argument)
output_dir="${@: -1}"

# Loop through each pattern directory
for motif_dir in "${motif_dirs[@]}"
do
    # Extract the motif name from the directory
    motif_name=$(basename "$motif_dir")

    # BED file for the current motif
    bed_files=(
        "${motif_dir}/msa/${motif_name}_CDS.bed"
        "${motif_dir}/msa/${motif_name}_nCDS.bed"
        "${motif_dir}/msa/${motif_name}_tIG.bed"
        "${motif_dir}/msa/${motif_name}_US.bed"
    )

  # Calling circular_plotter.R script
      Rscript "$circ" "${bed_files[0]}" "${bed_files[1]}" "${bed_files[2]}" "${bed_files[3]}"

done
