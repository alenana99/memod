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
> - R package ("GO.db")
> - R package ("optparse")
> - R package ("AnnotationDbi")
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

:warning: 	:warning:	 :warning: 	:warning:	 :warning: 	:warning: 	:warning: 	:warning: 	:warning:	 :warning: 	:warning: 	:warning:

# *This section is still a work in progress*

:construction: :construction: :construction: :construction: :construction: :construction: :construction: :construction: :construction: :construction: :construction: :construction: :construction: 

We are currently working on this part of the project. Check back soon for updates.

```
GSEA/
├── data/
│    ├── GO.out
│    ├── basecalling2_methylated_sites.tsv
├── scripts/
│    ├── process_meth.py
│    ├── process_go.R
│    ├── 
│    ├── me_GSEA.sh
└──  README.md
```

First, to do a GSEA, we need to get: 
- A gene list L = {g1,...gN} ranked according to the number of methylation (correlation r(g1)=r1)
- A gene set S (e.g. a pathway or a GO category)



```

```

```
```

```
```

```
```


```

```


```
res_ES <- lapply(list_S, function(x) ES_significance2(x,<L>))
```


```
NES_values <- sapply(myres_ES, function(x) x$NES)
NES_pi_values <- do.call(c, lapply(myres_ES, function(x) x$NES_pi))
```
## Reference
[*MeStudio* work](https://www.mdpi.com/1422-0067/24/1/159)

[*MicrobeMod* work](https://www.biorxiv.org/content/10.1101/2023.11.13.566931v1)

