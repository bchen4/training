
##ChIP-seq analysis pipeline

###Input files
* BAM file for each sample (From green center pipeline)
* BigWig file for each sample (From green center pipeline)
* Peak bed file for each sample (From green center pipeline)
* Design file
* Genome annotation files (including TSS, TTS, gene region bed files)

###Modules that can be processed parallel

####deeptTools
```bash
#Compute matrix using bigwig
computeMatrix reference-point -S bigwig_file -R Genome_range_bed_file -out matrix_file

#Draw region signal heatmap
plotHeatmap -m matrix_file --out heatmap.pdf

```
####DiffBind

```R
#Input is design file, R package will need bam files and peak bed files
R CMD run_diffbind.R design.txt

```
####ChIPSeeker

```R
R CMD run_chipseeker.R
```
