# MEMOD!
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white) ![Bash](https://img.shields.io/badge/bash-%234EAA25.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

Hi :blush:

Memod! is a workflow that uses multiple tools for exploring and analyzing prokaryotic methylation in Nanopore sequencing and combine genome-wide methylation profiles with genomic features. 


## Table of Contents

- [Dependencies](#everything-you-need-to-install-before-you-begin)
- [STEP 1: Basecalling with dorado](https://github.com/alenana99/memod/blob/main/README.md#step-1-basecalling-w-dorado)
- [STEP 2: Map reads to reference](https://github.com/alenana99/memod/blob/main/README.md#step-2-map-reads-to-reference)
- [STEP 3: MicrobeMod](https://github.com/alenana99/memod/blob/main/README.md#microbemod)
- [STEP 4: MeStudio](https://github.com/alenana99/memod/blob/main/README.md#mestudio)
- [STEP 5: Circular plots on R](https://github.com/alenana99/memod/blob/main/README.md#circular-plots-on-R)
- [STEP 6: Gene Set Enrichment Analysis](https://github.com/alenana99/memod/blob/main/README.md#gene-set-enrichment-analysis)
- [Reference](https://github.com/alenana99/memod/blob/main/README.md#reference)

# Everything you need to install before you begin
> - minimap2
> - samtools
> - R package ("optparse")
- **Dorado**
> - dorado
	
- **MicrobeMod**
 > - MicrobeMod
 > - Prodigal
 > - BLAST
 > - HMMER
 > - Modkit v0.2.2
 > - STREME
- **MeStudio**
 > - mestudio
 > - pandas
 > - numpy
 > - matplotlib
 > - R package ("circlize")
 > - readlink
 > - python 3.8 (or above)
 >
 - **GSEA**
 > - R package ("fgsea")
 > - R package ("dplyr")
 > - R package ("purrr")
## STEP 1: BASECALLING W/ DORADO

Dorado is a high-performance, easy-to-use, open source basecaller for Oxford Nanopore reads. 
See [Dorado](https://github.com/nanoporetech/dorado)

### Download methylation models
Check the [avaiable basecalling models](https://github.com/nanoporetech/dorado?tab=readme-ov-file#available-basecalling-models)
```
dorado download --model <CHOSEN MODEL>
```

### Convert fast5s in a pod5
If you have fast5, you can convert them in pod5 with the command like the below. See [Pod5 File Format Documentation](https://pod5-file-format.readthedocs.io/en/latest/docs/tools.html#pod5-convert-fast5)
```
pod5 convert fast5 ./*.fast5 --output converted.pod5
```

### Run dorado
```
dorado basecaller /dna_r10.4.1_e8.2_400bps_sup@v4.0.0 /converted.pod5 --modified-bases-models /dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v2,/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1 -x cpu > <LIBRARY_NAME>.bam 2> error_<LIBRARY_NAME>
```
This command uses a r10.4.1 basecalling model and passes two modified basecalling models as well, one for all context 4mC and 5mC and one for all context 6mA. The output is an unmapped BAM file with modified base information for each read.

## STEP 2: MAP READS TO REFERENCE

In order to map basecalled reads to the reference genome, including their methylation metadata, a bash-based script named *map_to_ref* has been implemented and you can use it. You can find the source code here.
*map_to_ref* expects three arguments:

1. <LIBRARY_NAME>.bam = BAM file obtained from Dorado basecalling
2.   reference.genomes.fna = the reference in FASTA format
3. <LIBRARY_NAME>.mapped.bam = output name

The flag -T (taglist) MM,ML in samtools fastq specify a comma separated list of methylation tags to copy to the FASTQ header line; minimap2 turns off secondary alignments; then, samtools view and sort generate a sorted bam and finally, the last samtools line is needed to index the BAM file.

## MicrobeMod

See [MicrobeMod](https://github.com/cultivarium/MicrobeMod) github repository for more informations

First, run MicrobeMod *call_methylation*
```
MicrobeMod call_methylation -b <LIBRARY_NAME>.mapped.bam -r reference_genomes.fna -t 10
```
We used the default parameters, which are robust in most cases, but it is possible to change them. You can re-running with different parameters to see which ones fit best with your data. Many files are created, of which the most important are the two tab-separated tables, one describing information for all methylated sites and one describing output for methylated motifs. These will be used later to continue the analysis with ```mestudio```.

Next, you can run MicrobeMod *annotate_rm*
```
MicrobeMod annotate_rm -f reference_genome.fna -o <reference_genome> -t 10
```
This pipeline just requires any genome assembly and creates many files, including a tabular output describing RM genes: each line describes an individual gene; genes are grouped in an "operon" (complete RM system) or a singleton (MTase or RE without a pair).

## MeStudio

MeStudio is a tool that enables the analysis and integration of genome-wide methylation profiles with genomic features. See the [MeStudio](https://github.com/combogenomics/MeStudio) github repository for more informations. We made just some adjustments to ensure compatibility with data from nanopore sequencing: we have modified the source codes for *ms_analyzR* and *imestudio*, while all the others remain unchanged. You can find the source codes here.

Before running ```imestudio```, make sure you have the files you need:
 - **<LIBRARY_NAME>_anno.gff** = the product of the genomic annotation, executed via ```prokka``` annotation; the file is reported in GFF3 format.
 - **<LIBRARY_NAME>_genomic.fasta** = the reference_genome.fna file, that you need to rename in order to make it usable for ```imestudio```
 - **motiflist.txt** = text file containing all the motifs you want to look for in the genome. You can use the list of motifs obtained from MicrobeMod *call_methylation* (in the output file <LIBRARY_NAME>_motif.tsv)
 - **<LIBRARY_NAME>_smart.gff** = GFF3 file with methylation positions. You can modify the MicrobeMod *call_methylation* output file <LIBRARY_NAME>.methylated_sites.tsv in order to obtain a GFF3 file with a similar format to that of the annotation file. 
 
 An example of how to obtain the format is provided below:

```
import pandas as pd
path_ms = "/<LIBRARY_NAME>_methylated_sites.tsv"
ms_df = pd.read_table(path_ms,sep="\t", low_memory=False)
ms_df = pd.DataFrame({
    'seqid': ms_df['Contig'],
    'source': 'dorado',
    'type': ms_df['Modification'],
    'start': ms_df['Position'],
    'end': ms_df['Position'],
    'score': '.',
    'strand': ms_df['Strand'],
    'phase': '.',
    'attribute': ms_df.apply(lambda row: f"coverage={row['Total_coverage']};context={row['Sequence']};Modified_bases={row['Modified_bases']};Unmodified_bases={row['Unmodified_bases']};Percent_modified={row['Percent_modified']};ID_gene={row['ID_gene']}", axis=1)
})
```

Manually add the path of the folder in which you want to locate ```imestudio``` by tiping: 
```
nano imestudio
```
In the INSTALLATION_DIR= (line 10) you can add your path and save your changes. If you don't have the administrator privileges, add the directory in which executables are present to the *$PATH*, for example: 
```
PATH=$PATH:/sources 
``` 

You have to render install and the other executables executable by typing: 
```
chmod +x <EXECUTABLE>
```

Now, you can run ```imestudio```
```
imestudio -f <LIBRARY_NAME_genomic.fasta> -anno <LIBRARY_NAME_anno.gff> -smart <LIBRARY_NAME_smart.gff> -mo <motiflist.txt> -o <output directory>
```

## Circular plots on R

For each searched motif, MeStudio generate a stdout log. ```mscore``` and ```msa``` are the directories created from ```ms_coreR``` and ```ms_analyzR```. Inside these folders you can find the tabular files as GFFs and BEDs produced by MeStudio. We have implemented the script *wrapper_motif.sh*, which internally calls the script *circular_plotter.R* to automatically generate circular plots for all searched motifs.

You can just type: 
```
bash wrapper_motif.sh <PATH/TO/MOTIF1> <PATH/TO/MOTIF2> <OUTPUT_DIRECTORY>
```
If you want, you can also match the BED differences: 
```
circos.initializeWithIdeogram(my_bed4)
colors <- c('#7fc97f','#beaed4','#fdc086','#ffff99')
bed_list <- list(my_bed, my_bed2, my_bed3, my_bed4)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = colors)
```
Feel free to change colors for the plot :wink:

![Our plot](https://github.com/alenana99/memod/blob/main/plots/plot_GATC_123.png)

## Gene Set Enrichment Analysis

```
GSEA/
├── data/
│    ├── GATC_CDS.bed
│    ├── GATC_nCDS.bed
│    ├── GATC_tIG.bed
│    ├── GATC_US.bed
│    ├── prokka_annotation
├── scripts/
│    ├── 
└──  README.md
```
Among the outputs generated from MeStudio, BED files were produced for each genomic feature (CDS, nCDS, tIG, and US) for each motif. These files contain columns representing the contig (seqid), start and end positions, the attribute (corresponding to the gene ID), the score (representing the number of methylations), and the gene annotation. The BED files are sorted in descending order by the score, with genes displaying the highest number of methylations listed first.


Let's merge MeStudio BED files output: GATC_CDS, GATC_nCDS, GATC_tIG, and GATC_US into one dataframe, replacing all NA values with 0:

```
GATC_total <- list(
  GATC_CDS %>% rename(score_CDS = score),
  GATC_nCDS %>% rename(score_nCDS = score),
  GATC_tIG %>% rename(score_tIG = score),
  GATC_US %>% rename(score_US = score)
) %>%
  reduce(function(x, y) merge(x, y, by = c("attribute", "seqid", "start", "end", "gene.annotation"), all = TRUE)) %>%
  mutate(across(everything(), ~ replace_na(., 0)))
```
To explore the relationship between gene length and the number of methylations, a LOESS regression model was applied. 

A new column 'r_exp' is created, containing the fitted values from the LOESS regression model, that represent the expected number of methylations based on gene length. A new column 'oe' (observed - expected) is created, which stores the difference between the observed and the expected number of methylations. This helps identify genes that have more or fewer methylations than expected based on their length. Another ew column 'obs/exp' is added, containing the ratio of observed to expected methylations. A value of 1 would indicate that the observed number of methylations matches the expected number, while values greater or less than 1 suggest over- or under-methylation relative to expected.

```
GATC_merged <- GATC_total %>%
  merge(ID_pos_length, by = "attribute") %>%
  mutate(
    r_exp = loess(score_CDS ~ length, data = .)$fitted,
    oe = score_CDS - r_exp,
    `obs/exp` = score_CDS / r_exp)
```

The following code creates a bar plot to visualize the residuals (i.e., the differences between observed and expected values), sorted in descending order. This allows you to quickly see which observations deviate the most from the LOESS model's predictions.
```
barplot(
  sort(GATC_merged$oe, decreasing = TRUE),
  main = "Plot of Sorted Residuals GATC",
  ylab = "Observed - Expected",
  border = "black",
  las = 2,  
  cex.names = 0.7  
)
```
Now, we are ready to perform GSEA with fgsea function from the FGSEA package. It requires: 
- pathways: A list of functional pathways (we use KEGG Pathways, but you can use others).
- stats: The ranked vector of observed-to-expected ratios.

To generate a list of pathways we can use the eggNOG-emapper annotation. This part creates a list of pathways where each KEGG pathway (KEGG_Pathway) is associated with the genes (query) that map to it. The split() function splits the genes vector by the KEGG_Pathway, producing a list where the names are the KEGG pathways and the elements are the genes belonging to each pathway.
```
pathways <- split(emapper.annotations$query, emapper.annotations$KEGG_Pathway)
```
To create the vector of gene statistics:
```
GATC_ranked_v <- setNames(GATC_tot$`obs/exp`, GATC_tot$attribute)
```
This creates a vector, where the values are the observed-to-expected ratios (obs/exp) from the dataset, and the names are gene identifiers (attribute). This vector will be used as the input for the GSEA analysis to rank the genes based on their methylation status.

Now, let's run fgsea:
```
GATC_fgseares <- fgsea(pathways = pathways,
		 stats = GATC_ranked_v,
		 scoreType = 'std',
		 nproc = 1)
```
scoreType = 'std': Indicates that both positive and negative scores are used in ranking.
nproc = 1: Runs the analysis using a single processor.

Extracting Significant Pathways
This section extracts the top 16 positively enriched pathways from the GSEA results. It filters pathways with positive NES (Normalized Enrichment Score), orders them by adjusted p-value (padj), and selects the top 16.

```
combined_pathways_GATC <- GATC_fgseares %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  head(16) %>%
  mutate(direction = "Positive") %>%
  bind_rows(
      GATC_fgseares %>%
      filter(NES < 0) %>%
      arrange(padj) %>%
      head(16) %>%
      mutate(direction = "Negative")
  ) %>%
  filter(pathway != "-")
```

This combines the positive and negative pathways into a single dataset and adds a new column direction to indicate whether the pathway is positively or negatively enriched.

To plot the NES values for significant pathways we can use ggplot. This generates a bar plot to visualize the NES values for the top 16 positively and negatively enriched pathways. The bars are colored based on the direction of enrichment:
- Positive pathways are filled with "skyblue".
- Negative pathways are filled with "deeppink".

The plot visually highlights which pathways are significantly enriched based on their NES values and adjusted p-values, making it easier to interpret which biological processes are over- or under-represented based on methylation levels.
```
ggplot(combined_pathways_GATC, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Positive" = "skyblue", "Negative" = "deeppink")) +
  coord_flip() +  
  labs(title = "Significant Pathways (Top 16)", x = "Pathway", y = "Adjusted P-value (padj)") +
  theme_minimal()
```


![Our plot](https://github.com/alenana99/memod/blob/main/plots/TopSignificantPathwaysGATC.pdf)

## Reference
[Riccardi, Christopher, Iacopo Passeri, Lisa Cangioli, Camilla Fagorzi, Marco Fondi, and Alessio Mengoni. 2023. "Crossing Bacterial Genomic Features and Methylation Patterns with MeStudio: An Epigenomic Analysis Tool" International Journal of Molecular Sciences 24, no. 1: 159. https://doi.org/10.3390/ijms24010159](https://www.mdpi.com/1422-0067/24/1/159)

[](https://www.biorxiv.org/content/10.1101/2023.11.13.566931v1)

