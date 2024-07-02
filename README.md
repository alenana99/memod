# MEMOD!
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
Hi :blush:
Memod! is a workflow that uses multiple tools for exploring and analyzing prokaryotic methylation in Nanopore sequencing and combine genome-wide methylation profiles with genomic features. 


## Table of Contents

- [Dependencies](#everything-you-need-to-install-before-you-begin)
- [STEP 1: Basecalling with dorado](https://github.com/alenana99/memod/blob/main/README.md#step-1-basecalling-w-dorado)
- STEP 2: Map reads to reference
- STEP 3: MicrobeMod
- STEP 4: MeStudio
- STEP 5: Circular plots on R
- STEP 6: Gene Set Enrichment Analysis
- Reference

# Everything you need to install before you begin
> - minimap2
> - samtools
> - R package ("GO.db")
- **Dorado**
> - dorado
	
- **MicrobeMod**
 > - MicrobeMod
 > - Prodigal
 > - BLAST
 > - HMMER
 > - Modkit v0.2.2
 > - STREME
- **mestudio**
 > - mestudio
 > - pandas
 > - numpy
 > - matplotlib
 > - R package ("circlize")
 > - readlink
 > - python 3.8 (or above)
 > 
## STEP 1: BASECALLING W/ DORADO

Dorado is a high-performance, easy-to-use, open source basecaller for Oxford Nanopore reads. 
See [Dorado](https://github.com/nanoporetech/dorado)

### Download methylation models
To check the [avaiable basecalling models](https://github.com/nanoporetech/dorado?tab=readme-ov-file#available-basecalling-models)
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
2.  reference.genomes.fna = the reference in FASTA format
3. <LIBRARY_NAME>.mapped.bam = output name

The flag -T (taglist) MM,ML in samtools fastq specify a comma separated list of methylation tags to copy to the FASTQ header line; minimap2 turns off secondary alignments; then, samtools view and sort generate a sorted bam and finally, the last samtools line is needed to index the BAM file.

## MicrobeMod

See [MicrobeMod](https://github.com/cultivarium/MicrobeMod) github repository for more informations

First, run MicrobeMod *call_methylation*
```
MicrobeMod call_methylation -b LIBRARY_NAME.mapped.bam -r reference_genomes.fna -t 10
```
I used the default parameters, which are robust in most cases, but it is possible to change them. You can re-running with different parameters to see which ones fit best with your data. Many files are created, of which the most important are the two tab-separated tables, one describing information for all methylated sites and one describing output for methylated motifs. these will be used later to continue the analysis with ```mestudio```.

Next, you can run MicrobeMod *annotate_rm*
```
MicrobeMod annotate_rm -f genome_reference.fasta -o genome_reference -t 10
```
This pipeline just requires any genome assembly and creates many files, including a tabular output describing RM genes: each line describes an individual gene; genes are grouped in an "operon" (complete RM system) or a singleton (MTase or RE without a pair).

## mestudio

Before running ```mestudio```, make sure you have the files you need:
 - **<LIBRARY_NAME>_anno.gff** = the product of the genomic annotation, executed via ```prokka``` annotation; the file is reported in GFF3 format.
 - **<LIBRARY_NAME>_genomic.fasta** = the reference_genome.fna file, that you need to rename in order to make it usable for ```mestudio```
 - **motiflist.txt** = text file containing all the motifs you want to look for in the genome. You can use the list of motifs obtained from MicrobeMod *call_methylation* (in the output file <LIBRARY_NAME>_motif.tsv)
 - **<LIBRARY_NAME>_smart.gff** = GFF3 file with methylation positions. You can modify the MicrobeMod *call_methylation* output file <LIBRARY_NAME>.methylated_sites.tsv in order to obtain a GFF3 file with a similar format to that of the annotation file. 
 
 An example of the format is reported below:

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

Manually add the path of the folder in which you want to locate ```mestudio``` by tiping: 
```
nano imestudio
```
In the INSTALLATION_DIR= (line 10) you can add your path and save your changes. If you don't have the administrator privileges, add the directory in which executables are present to the *$PATH*, for example: 
```
PATH=$PATH:/src 
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

Import library
```
library(circlize)
```

Read your bed files (imestudio output files)
```
my_bed <- as.data.frame(read.table("[MOTIF]_CDS.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
my_bed2 <- as.data.frame(read.table("[MOTIF]_nCDS.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
my_bed3 <- as.data.frame(read.table("[MOTIF]_intergenic.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
my_bed4 <- as.data.frame(read.table("[MOTIF]_upstream.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
```
Init plot:
```
circos.initializeWithIdeogram(my_bed4)
```
Specify each bed's methylation-density: 
```
circos.genomicDensity(my_bed, col = c("red"), track.height = 0.1)
circos.genomicDensity(my_bed2, col = c("blue"), track.height = 0.1)
circos.genomicDensity(my_bed3, col = c("purple"), track.height = 0.1)
circos.genomicDensity(my_bed4, col = c("orange"), track.height = 0.1)
```
If you want to match the bed differences:
```
colors <- c('#7fc97f','#beaed4','#fdc086','#ffff99')
bed_list = list(my_bed, my_bed2, my_bed3, my_bed4)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = colors)
```
## Gene Set Enrichment Analysis

First, to do a GSEA, we need to get: 
- A gene list L = {g1,...gN} ranked according to the number of methylation (correlation r(g1)=r1)
- A gene set S (e.g. a pathway or a GO category)

We want to associate a GO term with each gene; to do this, we can use [PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/), an automatic service for functional annotation of prokaryotic and eukaryotic proteins of unknown function. We upload the protein file obtained from prokka annotation to PANNZER. In the parameters filters we can select Bacteria to exclude unlikely GO terms. PANNZER2 return us a file <GO.out> with all the GO terms associated with a certain gene of our gene universe, that is, of all the genes that we take into account. 

To get a list of all GO ids:
```
cut -f3 GO.out > GO_ids
```
Now, we can use *fetch_go_details* function of R in order to obtain a dataframe with ontology, name and definition for each GO term: 
```
library(GO.db)
GO_ids <- read.table("GO_ids")
GO_ids <- as.vector(as.matrix(GO_ids)
GO_description <- fetch_go_details(GO_ids)
```

Let's merge each of our gene ids with its GO informations: 
```
cut -f1,3 GO.out > geneID_goID
geneID_goID <- read.csv("~/myproggg/test/geneID_goID", sep="")
names(geneID_goID) <- c("gene_ID", "GO_ID")
geneID_goID_merged <- merge(geneID_goID, GO_description, by = "GO_ID")
geneID_goID_merged <- unique(geneID_goID_merged)
```
Let's divide by categories: Molecular Function, Cellular Component and Biological Process and then construct an S for each GO of each category:
```
BP <- subset(geneID_goID_merged, Ontology == "BP")
MF <- subset(geneID_goID_merged, Ontology == "MF")
CC <- subset(geneID_goID_merged, Ontology == "CC")
myS_BP <- data.frame(BP$gene_ID)
names(myS_BP) <- c("g")
```
Now, we can use the methylated positions from the MicrobeMod call_methylation output file (<LIBRARY_NAME>_methylated_sites.tsv) and the start and end positions of each gene from the GFF3 annotation file to count the number of methylations per gene.
!! We can also first divide by methylation type and then count the number of methylations and, if we want, we can do the same for upstream positions (200bp upstream)
```
import pandas as pd
path_pos_ID = "/pos_ID"
pos_ID_df = pd.read.table(path_pos_ID, header = None)
pos_ID_df.columns = ['pos_in', 'pos_end', 'ID_gene']
path_pmsm = "/pos_mod_strand_motif.tsv"
pmsm_df = pd.read.table(path_pmsm)

def find_gene(position, pos_ID_df):
    for _, row in pos_ID_df.iterrows():
        if row['pos_in'] <= position <= row['pos_end']:
            return row['ID_gene']
    return None

pmsm_df['ID_gene'] = pmsm_df['Position'].apply(lambda pos: find_gene(pos,pos_ID_df))
pmsm_df['ID_gene'] = pmsm_df['ID_gene'].str.replace('_gene','')

pmsm_count_df = pmsm_df.groupby(['ID_gene'])['ID_gene'].count()
```
Sort by number of methylations:
```
pmsm_count_df = pmsm_count_df.sort_values(ascending=False)
pmsm_count_df = pmsm_count_df.reset_index(name='n_meth')
pmsm_count_df['ID_gene'] = pmsm_count_df['ID_gene'].str.replace('ID=', '')
```
```
pmsm_count_df <- read.table("/pmsm_count_df", header=TRUE)
S_BP$i <- match(S_BP$g, pmsm_count_df$g)
S_BP_GO <- data.frame(BP$gene_ID, BP$GO_ID)
names(S_BP_GO) <- c("g", "GO_ID")
S_BP_GO$i <- match(S_BP_GO$g, pmsm_count_df$g)
S_BP_GO <- merge(S_BP_GO, pmsm_count_df, by = "g")
unique_GO_BP <- unique(S_BP_GO$GO_ID)
```
We can create a list containing $g (gene_ID) and $i ( match S\$g, L$g)  for each GO term of each category (now we use BP for example):
```
list_S <- list()
for (go_id in myunique_GO_BP) {
  ss <- myS_BP_GO[myS_BP_GO$GO_ID == go_id, c("g", "i", "r")]
  ss <- ss[order(ss$r, decreasing = TRUE), ]
  ss <- ss[order(ss$i), ]
  ss <- ss[, c("g", "i")]
  list_S[[paste0("", gsub(":", "_", go_id))]] <- ss
}
```
Finally, we are ready for GSEA!
First, with *ES_significance2* function we calculate an Enrichment Score (ES, with *ES_wrapper_fast*) that reflect the degree to which a set S is overrapresented at the top or bottom of the entire ranked list L. The score is calculated by waling down the list L, increasing a running-sum statistic (Phit) when we encounter a gene in S and decreasing it when we encounter genes not in S (on the contrary, Pmiss increases). Then it estimates the statistical significance (P value) of the ES by using a permutation test (bootstrap = 1000): at each iteration, we randomly sample dim(x)[1] rows from L and calculate the ES for each sample and store the result in *ES_pi*. After calculating significance, we adjust the estimated significance level to account for multiple hypothesis testing. We first normalize the ES for each gene set to account for the size of the set, yielding a normalized enrichment score (NES).
The function return a list containing for each set:
- **ES**: the initially calculated enrichment statistic
- **sig**: the calculated significance
- **NES**: the normalized enrichment statistic
- **NES_pi**: the normalized bootstrap values

```
res_ES <- lapply(list_S, function(x) ES_significance2(x,pmsm_count_df))
```

From the obtained list we extract **NES** and **NES_pi** and use them as inputs for *FDR_wrapper* function to compute the False Discovery Rates (FDRs). The FDR is the estimated probability that a set with a given NES rapresents a false positive finding. 
```
NES_values <- sapply(myres_ES, function(x) x$NES)
NES_pi_values <- do.call(c, lapply(myres_ES, function(x) x$NES_pi))
```
