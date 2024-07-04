library(circlize)

# Ottieni i file BED dai parametri della riga di comando
args <- commandArgs(TRUE)

# Lettura dei file BED
my_bed <- as.data.frame(read.table(args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed2 <- as.data.frame(read.table(args[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed3 <- as.data.frame(read.table(args[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed4 <- as.data.frame(read.table(args[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))

# Inizializzazione del plot circolare
circos.initializeWithIdeogram(my_bed4)

# Specifica della densit▒|  di metilazione per ciascun BED
circos.genomicDensity(my_bed, col = c("#7fc97f"), track.height = 0.1)
circos.genomicDensity(my_bed2, col = c("#beaed4"), track.height = 0.1)
circos.genomicDensity(my_bed3, col = c("#fdc086"), track.height = 0.1)
circos.genomicDensity(my_bed4, col = c("#ffff99"), track.height = 0.1)

# Se vuoi abbinare le differenze BED
#colors <- c('#7fc97f','#beaed4','#fdc086','#ffff99')
#bed_list <- list(my_bed, my_bed2, my_bed3, my_bed4)
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = colors)
