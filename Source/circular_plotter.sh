#!/bin/bash

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

    # Cambia la directory alla cartella "msa" all'interno del motivo
    cd "$motif_dir/msa" || { echo "Error: cannot change directory to $motif_dir/msa"; exit 1; }

    # File BED per il motivo corrente
    bed_files=(
        "${motif_name}_CDS.bed"
        "${motif_name}_nCDS.bed"
        "${motif_name}_intergenic.bed"
        "${motif_name}_upstream.bed"
    )

    # Creazione dello script R per i plot circolari
    Rscript - <<EOF
library(circlize)

# Lettura dei file BED
my_bed <- as.data.frame(read.table("$bed_files[1]", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed2 <- as.data.frame(read.table("$bed_files[2]", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed3 <- as.data.frame(read.table("$bed_files[3]", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed4 <- as.data.frame(read.table("$bed_files[4]", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))

# Inizializzazione del plot circolare
circos.initializeWithIdeogram(my_bed4)

# Specifica della densità di metilazione per ciascun BED
circos.genomicDensity(my_bed, col = c("red"), track.height = 0.1)
circos.genomicDensity(my_bed2, col = c("blue"), track.height = 0.1)
circos.genomicDensity(my_bed3, col = c("purple"), track.height = 0.1)
circos.genomicDensity(my_bed4, col = c("orange"), track.height = 0.1)

# Se vuoi abbinare le differenze BED
colors <- c('#7fc97f','#beaed4','#fdc086','#ffff99')
bed_list <- list(my_bed, my_bed2, my_bed3, my_bed4)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = colors)

# Salva il plot
pdf("$output_dir/$motif_name\_circos_plot.pdf")
circos.clear()
EOF

    echo "Plot circolare per il motivo $motif_name creato: $output_dir/$motif_name\_circos_plot.pdf"

    # Torna alla directory di partenza dopo aver completato il lavoro su questo motivo
    cd - >/dev/null || exit
done
