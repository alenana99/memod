library(circlize)

# Get BED files from command line parameters
args <- commandArgs(TRUE)

# Reading BED files
my_bed <- as.data.frame(read.table(args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed2 <- as.data.frame(read.table(args[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed3 <- as.data.frame(read.table(args[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
my_bed4 <- as.data.frame(read.table(args[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))

# Initialization of the circular plot
circos.initializeWithIdeogram(my_bed4)

# Specification of methylation density for each BED
circos.genomicDensity(my_bed, col = c("#7fc97f"), track.height = 0.1)
circos.genomicDensity(my_bed2, col = c("#beaed4"), track.height = 0.1)
circos.genomicDensity(my_bed3, col = c("#FFFF7F50"), track.height = 0.1)
circos.genomicDensity(my_bed3, col = c("#FF6495ED"), track.height = 0.1)

# If you want to match BED differences
#circos.initializeWithIdeogram(my_bed4)
#colors <- c('#7fc97f','#beaed4','#fdc086','#ffff99')
#bed_list <- list(my_bed, my_bed2, my_bed3, my_bed4)
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = colors)
