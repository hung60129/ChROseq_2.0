# ChROseq_2.0

Instructions on running ChRO-seq jobs in [the Sethupathy Lab](https://github.com/Sethupathy-Lab). See [getting ready to run a job](https://github.com/Sethupathy-Lab/cornell_tutorials/blob/master/getting_ready_to_run_a_job.md) for general guidance of running josb at Cornell bioHPC environment. <br>

For more in-depth information about ChRO-seq method development, check out [Chu et al., Nature Genetics (2018)](https://www.nature.com/articles/s41588-018-0244-3) published by the Danko lab. A few more examples of how ChRO-seq can be applied to study human health and diseases: [Dinh et al., Cell Reports (2020)](https://www.sciencedirect.com/science/article/pii/S2211124720303995) for studying a rare liver cancer type, [Hung et al., NAR (2021)](https://academic.oup.com/nar/article/49/2/726/6066631) for studying human gut development, [Hung et al., bioRxiv](https://www.biorxiv.org/content/10.1101/2021.08.20.457162v2) for studying non-alcoholic steatohepatitis (NASH), and [Hung et al., bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.12.499825v1.full) for a multi-omics integration study. 

## 1. Genome mapping
>**Goal**: Eliminate PCR duplicates, trim reads, QC, map to a copy of the genome that contains the rRNA PolI repeating transcriptional unit. <br>
>**Tool source**: Forked from https://github.com/Danko-Lab/proseq2.0/blob/master/proseq2.0.bsh with minor modifications. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper` <br>
>**Appropriate compute node**: 40 gpu

Prepare mapping files (below are the path of genomes in the Sethupathy lab): <br>
   - mouse (mm9): `/home/pr46_0001/projects/genome/mm9_rRNA` <br>
   - human (hg38): `/home/pr46_0001/projects/genome/GRCh38.p7_rRNA` <br>
   - rat (rn6): `/home/pr46_0001/projects/genome/rn6_rRNA` <br>
   - If working with other species, use [bwa command](http://bio-bwa.sourceforge.net/bwa.shtml) to generate mapping files, add rRNA sequence to `.fa` file, and modify chromosom size file accordingly. <br>

Once on a compute node, source the environment to run ChRO-seq mapping: <br>
```
source /home/pr46_0001/ChROseq/ChROseq_pipeline/ProseqMapper/setChROenv.bsh
```

In the directory containing fastq.gz files, execute the job: <br>
```
bash /home/pr46_0001/ChROseq/ChROseq_pipeline/ProseqMapper/proseq2.0.bsh \
  -SE \
  -G \
  -i /home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.primary_assembly.genome_rRNA \
  -c /home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.rRNA.chrom.sizes \
  -O [PROJECT_NAME]_ProseqMapper2.0_output \
  --UMI1=6 \
  --thread=32 \
  -4DREG &> Proseq_output.log&
```

In the output folder, you should find 5 output files per sample, `proseq2.0_read_report.log`, and `proseq2.0_Run.log`. <br>
```
[sample_A]_dedup_QC.align.log
[sample_A].prinseq-pcrDups.gd
[sample_A]_dedup_QC_minus.bw
[sample_A]_dedup_QC_plus.bw 
[sample_A]_dedup_QC.sort.bam
proseq2.0_read_report_XXXXXXXXXXXXXXXXXXXXXXXXXX.log
proseq2.0_Run_XXXXXXXXXXXXXXXXXXXXXXXXXXX.log
```

You will collect mapping statistics in this output folder: <br>
```
python3 /workdir/your_cornell_ID/ProseqMapper/collect_ChROseq_mapping_statistics.py \
  [PROJECT_NAME]_mapping_statistics.txt \
  proseq2.0_read_report_XXXXXXXXXXXXXXXXXXXXXXXXXXXX.log
```

*Notes*: % Reads mapped should be ~5-35% (at least 20% preferred).

## 2. Merging bigwig files
>**Goal**: Merge bigwig files from mapping for the next step of TRE calling  <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper/mergeBigWigs.bsh` <br>
>**Appropriate compute node**: 40 gpu <br>

Set path to `kentUtils` tool: <br>
```
export PATH=/programs/kentUtils/bin:$PATH
```

In mapping output folder, merge mapping from all plus.bw files: <br>
```
bash /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper/mergeBigWigs.bsh \
  -c /home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.rRNA.chrom.sizes [PROJECT_NAME]_ALL_plus.bw *_plus.bw
```

Repeat the step to merge minus.bw files: <br>
```
bash /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper/mergeBigWigs.bsh \
  -c /home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.rRNA.chrom.sizes [PROJECT_NAME]_ALL_minus.bw *_minus.bw
```

## 3. TRE calling
>**Goal**: Define a universal set of transcriptional regulatory elements (TREs) using the latest [dREG tool](https://github.com/Danko-Lab/dREG/) from the Danko lab. Essentially, dREG is a machine learning algorithm that search for short bi-directional expression pattern of ChRO-seq signal across the entire genome. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/dREG/run_dREG_multiSample.bsh` <br>
>**Appropriate compute node**: 40 gpu (*essential to use gpu for this machine learning algorithm!*) <br>

Edit bash script for `outdir` and `gpu` parameters (if needed): <br>
*Note*: do not change `asvm` parameter. That file is the original dataset used to train dREG. <br>
```
nano run_dREG_multisample.bsh

########################################
###### Input dREG parameters here ######
########################################

outdir='dREG_calling_output'
asvm='/home/pr46_0001/ChROseq/ChROseq_pipeline/dREG/dREG_models/run_dREG/asvm.gdm.6.6M.20170828.rdata'
cores=30 #number of cpu cores
gpu=0 #use 0 if running with gpu, use '' if running without gpu

########################################
```

Execute the job: <br>
```
bash /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/dREG/run_dREG_multiSample.bsh ./[PROJECT_NAME]_mapping_output/*_plus.bw
```

The dREG output directory will show up after the job is complete. This command generates bed files with peak annotation, scores, and additional information for each individual sample and the “merged” bigwig files. Detailed output information is available [here](https://github.com/Danko-Lab/dREG/). We used `*.dREG.peak.score.bed.gz` files for downstream analysis. <br>

In dREG output directory, you can do a quick check for # of TREs identified in your samples: <br>                                                         
```
for i in *peak.score.bed.gz; do echo $i; cat $i | gunzip | wc -l; done 
```

After this step, we switch from a GPU to a standard 24-core compute node. <br>


## 4. Data visualization
### Bigwig file normalization
>**Goal**: Perform DESeq2-based differential expression analysis of TREs. <br>
>**R script template**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ChROseq_R/DESeq2_ChROseq_intragenicTRE_tutorial.R` <br>
>**Appropriate compute node**: 24-core node <br>

### Greate a genome track on UCSC genome browser
Follow 
To generate genome tracks for your project, 

### Other tools
You can consider using tools such as `DeepTool` to visualize signal distribution. See an example from [my study]() - Supplementary figure 1. 


## 5. TRE type classification
>**Goal**: Sort TREs into enhancer or promoter based on their genomic coordinates using a in-house tool. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools_intragenicTRE/classifyTRE_intragenicTRE_v2.0.py` <br>
>**Appropriate compute node**: 24-core node <br>

You will need a gtf file for gene coordinate information (below are the path of genomes in the Sethupathy lab): <br>
- mouse (mm9): `/home/pr46_0001/projects/genome/mm9/gencode.vM1.annotation.gtf` <br>
- human (hg38): `/home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf` <br>
- rat (rn6): `/home/pr46_0001/projects/genome/rn6_rRNA/Rattus_norvegicus.Rnor_6.0.91.gtf` <br>

Subject `dREG.peak.score.bed.gz` file to the in-house python tool. The example command defined window of 1000 base upstream and 200 base downstream of TSS as a promoter region.  <br>
```
zcat [PROJECT_NAME]_ALL.dREG.peak.score.bed.gz | python3 classifyTRE_intragenicTRE_v2.0.py \
 -u 1000 \ # upstream window
 -d 200 \ # downstream window
 -o 1 \ # number of base to be considered 'overlapping'
 /home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf \
 > [PROJECT_NAME]_classifyTRE_intragenicTRE_v2.0_output.txt 
```

How to inteprete output: TREs may be part of multiple genes (e.g. overlapping genes), multiple promoters (e.g. bidirectional promoters), or both promoter and genic (e.g. immediate downtream of a transcription start site). All possibilities are output here. Output includes the TRE bed file followed by the columns below. <br>
- **genic status**: whether a query TRE overlaps with any annotated gene bodies (TRUE or FALSE) <br>
- **overlapping genes**: If genic status == TRUE, what is the gene name(s)? <br>
- **proximal status**: whether a query TRE overlaps with any promoter windows (TRUE or FALSE) <br>
- **overlapping transcript TSSs**: If proximal status == TRUE, what is the gene name(s) of the TSS? <br>
- **distal status**: TRUE if a query TRE is not overlapped with any gene or transcript coordinates; considered as an ehancer. <br>
- **intragenic status**: whether a query TRE is entirely within a gene body (TRUE or FALSE) <br>
- **strand**: If intragenic status == TRUE, which strand of the overlapped gene body (+ or -) <br>

Do a quick check for numbers of enhancers and promoters that are detected in your sample set. <br>
```
awk '$9 == "Distal:True" {print $0}' HP_classifyTRE_v2.0_batch12_TREs.txt | wc -l
awk '$9 == "Distal:FALSE" {print $0}' HP_classifyTRE_v2.0_batch12_TREs.txt | wc -l
```
You should expect to see enhancer #: promoter # ratio close to 2:1. However, this may depends on cell/tissue types and mapping rates. The number of detectable enhancer TREs is decreased if the sample has low mapping coverage.   

*Alternatives*: There are other genomic tools out there that can annotate TRE types. Below is an example of using [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html) to run this job.
```
source /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/HOMER/setHOMERenv.bsh

annotatePeaks.pl [PROJECT_NAME]_ALL.dREG.peak.score.bed hg38 > [PROJECT_NAME]_HOMER_TRE_annotation_output.txt
```

## 6. Define differentially transcribed TREs between cell types/conditions
>**Goal**: Perform DESeq2-based differential expression analysis of TREs. <br>
>**R script template**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ChROseq_R/DESeq2_ChROseq_intragenicTRE_tutorial.R` <br>
>**Appropriate compute node**: 24-core node <br>

See detailed instructions [here](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=266#c) for using Rstudio at Cornell BioHPC. Make sure to load R 3.5.0 for this job. <br>
```
/programs/rstudio_server/rstudio_start 3.5.0
```

You can modify the R template for your project needs. In this latest version of script, we quantify TRE signal differently based on the `intragenic status` of TRE type output. **For those that are entirely within gene bodies (intragenic TREs), we extract counts only from the opposite (non-stranded) strand and multiply the singal by 2. For those that are not entirely within gene bodies (non-intragenic TREs), we extract and sum counts from both strands.** This change was made for eliminating bias from gene transcription activity. <br>


## 7. Motif enrichment analysis


## 8. Super-enhancer analysis


## 9. Define differentially transcribed genes between cell types/conditions


## 10. Other ChRO-seq related tools
### Define de novo transcription units

### Histome modificaiton calling


