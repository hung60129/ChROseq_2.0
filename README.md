# ChROseq_2.0

Instructions on running ChRO-seq jobs in [the Sethupathy Lab](https://github.com/Sethupathy-Lab). See [getting ready to run a job](https://github.com/Sethupathy-Lab/cornell_tutorials/blob/master/getting_ready_to_run_a_job.md) for general guidance of running josb at Cornell bioHPC environment. For more in-depth information about ChRO-seq method development, check out [Chu et al., Nature Genetics (2018)](https://www.nature.com/articles/s41588-018-0244-3) published by the Danko lab. 

ChROseq analysis v1.0 was built by the previous member Dr. Tim Dihn. ChRO-seq v1.0 was described in [Dinh et al., Cell Reports (2020)](https://www.sciencedirect.com/science/article/pii/S2211124720303995) and [Hung et al., NAR (2021)](https://academic.oup.com/nar/article/49/2/726/6066631). I have made several changes to make this ChROseq_2.0 pipeline and applied it to studies [Hung et al., bioRxiv (NASH project)](https://www.biorxiv.org/content/10.1101/2021.08.20.457162v2) and [Hung et al., bioRxiv (multi-omics project)](https://www.biorxiv.org/content/10.1101/2022.07.12.499825v1.full). 


## 1. Genome mapping
>**Goal**: Eliminate PCR duplicates, trim reads, QC, map to a copy of the genome that contains the rRNA PolI repeating transcriptional unit. <br>
>**Tool source**: Forked from https://github.com/Danko-Lab/proseq2.0/blob/master/proseq2.0.bsh with minor modifications. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper` <br>
>**Appropriate type of compute node**: 40 gpu

Prepare mapping files (below are the path of genomes in the Sethupathy lab): <br>
   - mouse (mm9): `/home/pr46_0001/projects/genome/mm9_rRNA` <br>
   - human (hg38): `/home/pr46_0001/projects/genome/GRCh38.p7_rRNA` <br>
   - rat (rn6): `/home/pr46_0001/projects/genome/rn6_rRNA` <br>

If working with other species, use [bwa command](http://bio-bwa.sourceforge.net/bwa.shtml) to generate mapping files, add rRNA sequence to `.fa` file, and modify chromosom size file accordingly. Also, currently we only has mm9 genome build for mouse samples. You can consider use `liftOver` to convert mm9 output to mm10 if nedded. <br>

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


## 2. TRE calling
### 2-1: Merge bigwig files 
>**Goal**: Merge bigwig files from mapping for the next step of TRE calling  <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper/mergeBigWigs.bsh` <br>
>**Appropriate type of compute node**: 40 gpu <br>

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

### 2-2: TRE calling using dREG
>**Goal**: Define a universal set of transcriptional regulatory elements (TREs) using the latest [dREG tool](https://github.com/Danko-Lab/dREG/) from the Danko lab. Essentially, dREG is a machine learning algorithm that search for short bi-directional expression pattern of ChRO-seq signal across the entire genome. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/dREG/run_dREG_multiSample.bsh` <br>
>**Appropriate type of compute node**: 40 gpu (*essential to use gpu for this machine learning algorithm!*) <br>

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


## 3. Data visualization
### 3-1: Merge Bigwig files for visualization purpose
>**Goal**: Merge plus & minus bigwig files from mapping step by experimental group for visualization on the UCSC genome browser as individual tracks. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/normalize_stranded_bigwigs4visualization.sh` <br>
>**Appropriate type of compute node**: 24-core node <br>

Noted that this bigwig merging job is different from the one in the TRE calling step. For TRE calling step, we merge bigwig files from ALL the samples WITHOUT normalization. For visualization purpose, we merge bigwig files for each of the experimental groups followed by a normalization process. <br>

Some notes from Tim: "Bigwig signal is determined by two things: (1) number of reads (2) read length. In some cases, read length is variable (smRNA-seq, ChRO-seq). In these situations, normalization is tricky because you can't normalize accurate using a per read metric (unless you get an average read length like RNA-seq or represent each read as a single value regardless of read length). Here we normalize using total bigwig signal, "wigsum". This script normalizes to a total wigsum of 100,000,000." <br>

First, source the shell script and tool path by the following commands: <br>
```
source /programs/RSeQC2-2.6.1/setup.sh
export PATH=/programs/kentUtils/bin:$PATH
```

Next, edit the bash script related to chromosome path (the example below is for human samples): <br>
```
nano normalize_stranded_bigwigs4visualization.sh

####################################
###### Input  parameters here ######
####################################

chrom='/home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.rRNA.chrom.sizes'
wigsum=100000000
suffix='normalized_100Msignal'
```
Also further down in the script: <br>
```
pos_sum=`cat tmp_plus.bedGraph | sort -k1,1 -k2,2n | bedtools map -a /home/pr46_0001/projects/genome/GRCh38.p7/GRCh38.chrom.sizes.bed -b stdin -c 4 -o sum | awk '{sum+=$4} END {print sum}'`
neg_sum=`cat tmp_minus.bedGraph | sort -k1,1 -k2,2n | bedtools map -a /home/pr46_0001/projects/genome/GRCh38.p7/GRCh38.chrom.sizes.bed -b stdin -c 4 -o sum | awk '{sum+=$4} END {print sum}'`
```

Finally, to carry out this merging+normalization step, you merge bigwigs files using `mergeBigWigs.bsh` followed by `normalize_stranded_bigwigs4visualization.sh`. Be noted that you tackle plus and minus bigwig files saparately and tackle experimental groups separately. <br>
```
bash /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper/mergeBigWigs.bsh \
  -c /home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.rRNA.chrom.sizes [PROJECT_NAME]_combined_StageA_plus.bw StageA_sample1_plus.bw StageA_sample2_plus.bw StageA_sample3_plus.bw
```
```
bash /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ProseqMapper/mergeBigWigs.bsh \
  -c /home/pr46_0001/projects/genome/GRCh38.p7_rRNA/GRCh38.rRNA.chrom.sizes [PROJECT_NAME]_combined_StageA_minus.bw StageA_sample1_minus.bw StageA_sample2_minus.bw StageA_sample3_minus.bw
```
You input the path to `[PROJECT_NAME]_combined_StageA_plus.bw` when executing the normalization bash file. It will automatically do the normalization for the minus bigwig file if it has the same file prefix. <br>
```
bash /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/normalize_stranded_bigwigs4visualization.sh \  
[PROJECT_NAME]_combined_StageA_plus.bw &> bigwig_visualization_output.log&
```
The complation of the job will generate files with suffix 'normalized_100Msignal'. You will upload the normalized bigwig files to UCSC genome browser. <br>

### 3-2: Greate a genome track on UCSC genome browser
You can create tracks to visuzlize TREs and ChRO-seq intensity track on the plus/minus strand to [UCSC genome browser](https://genome.ucsc.edu/). We use ftp to upload files. Find more instructions on UCSC genome browser for file uplad and track formatting. Below are the general guideline you can follow:   

1. Log into UCSC account and upload files as “Custom Tracks” <br>
2. Choose appropriate genome assembly: Human (hg38) → “Add Custom Track” <br>
3. Add and submit URLs to files individually (My Data → Custom Tracks). Below are some examples: <br>
```
ftp://cbsuftp.biohpc.cornell.edu/pr46ftp/[YOUR_CORNELL_ID]/[PROJECT_NAME]_all.dREG.peak.score.bed
ftp://cbsuftp.biohpc.cornell.edu/pr46ftp/[YOUR_CORNELL_ID]/[PROJECT_NAME]_combined_stageA_plus_normalized_100Msignal.bw 
ftp://cbsuftp.biohpc.cornell.edu/pr46ftp/[YOUR_CORNELL_ID]/[PROJECT_NAME]_combined_stageA_minus_normalized_100Msignal.bw
```
4. View in genome browser → “Go” <br>
5. Check some genes of interest to make sure files are correct <br>
6. Save Settings to allow session to be loaded and shared with others. “My Session” → “Save Settings” <br>
7. You can configure settings to edit layout – any updates need to be saved <br>
8. You can share URL to collaborators! <br>

I like to add the following line at the very beginning of bed files of TRE coordinates before uploaded onto UCSC Genome Browser: add custom track name and description if you like. <br>
```
track name=dREG_peak_score_ALL_TREs description="dREG peak score (ALL TREs)" itemRgb=On
```

### 3-3: Other visuzalization tools 
You can consider using tools such as [DeepTool](https://deeptools.readthedocs.io/en/develop/) to visualize signal distribution. See an example from [here](https://www.biorxiv.org/content/10.1101/2022.07.12.499825v1.full) - Supplementary figure 1. <br>


## 4. TRE type classification
>**Goal**: Sort TREs into enhancer or promoter based on their genomic coordinates using a in-house tool. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools_intragenicTRE/classifyTRE_intragenicTRE_v2.0.py` <br>
>**Appropriate type of compute node**: 24-core node <br>

You will need a gtf file for gene coordinate information (below are the path of genomes in the Sethupathy lab): <br>
- mouse (mm9): `/home/pr46_0001/projects/genome/mm9/gencode.vM1.annotation.gtf` <br>
- human (hg38): `/home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf` <br>
- rat (rn6): `/home/pr46_0001/projects/genome/rn6_rRNA/Rattus_norvegicus.Rnor_6.0.91.gtf` <br>

Subject `dREG.peak.score.bed.gz` file to the in-house python tool. The example command defined window of 1000 base upstream and 200 base downstream of TSS as a promoter region.  <br>
```
zcat [PROJECT_NAME]_ALL.dREG.peak.score.bed.gz | python3 classifyTRE_intragenicTRE_v2.0.py \
 -u 1000 \
 -d 200 \
 -o 1 \ 
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


## 5. Define differentially transcribed TREs between cell types/conditions
>**Goal**: Perform DESeq2-based differential expression analysis of TREs. <br>
>**R script template**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ChROseq_R/DESeq2_ChROseq_intragenicTRE_tutorial.R` <br>
>**Appropriate type of compute node**: 24-core node <br>

See detailed instructions [here](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=266#c) for using Rstudio at Cornell BioHPC. Make sure to load R 3.5.0 for this job. <br>
```
/programs/rstudio_server/rstudio_start 3.5.0
```

You can modify the R template for your project needs. The R script include the 3 major steps (see below). <br>
1. Extract ChRO-seq counts from TRE regions using package `bigwig` <br>
2. Exploratory data analysis including heirachical clustering analysis and data dimentionality reduction by PCA 
3. Input count matrix for differential expression analysis using package `DESeq2` <br>

In this latest version of script, we quantify TRE signal differently based on the `intragenic status` of TRE type output. **For those that are entirely within gene bodies (intragenic TREs), we extract counts only from the opposite (non-stranded) strand and multiply the singal by 2. For those that are not entirely within gene bodies (non-intragenic TREs), we extract and sum counts from both strands.** This change was made for eliminating bias from gene transcription activity. <br>


## 6. Motif enrichment analysis
>**Goal**: Use tool [HOMER](http://homer.ucsd.edu/homer/index.html) to identify transcription factor binding motifs enriched in a given set of TRE regions. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/HOMER/` <br>
>**Appropriate type of compute node**: 24-core node <br>

### 6-1: Prepare bed files
You can use the following R template to generate files of TREs of interest and proper background regions in bed format. <br>  
`/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ChROseq_R/HOMER_analysis_intragenicTRE_tutorial.R` <br>

See detailed instructions [here](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=266#c) for using Rstudio at Cornell BioHPC. Make sure to load R 3.5.0 for this job. <br>
```
/programs/rstudio_server/rstudio_start 3.5.0
```

Say here you want to compare TREs that are specific to stage A compared to control. You can use this template to generate `[PROJECT-NAME-stageA-specific-TRE].bed` and `[PROJECT-NAME-non-stageA-specific-TRE].bed` (background), both will be used for HOMER motif analysis. <br>

### 6-2: HOMER motif analysis
Find more information about HOMER motif analysis [here](http://homer.ucsd.edu/homer/motif/). This tool has been installed in our lab space, all you need to do is to source the HOMER environment by the following command and you should be able to carry out any analyses provided by HOMER.  <br>
```
source /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/HOMER/setHOMERenv.bsh
```

We use function `findMotifsGenome.pl` from HOMER. Usage of the command is: `findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]`. We execute this job (see the example command below) by providing TREs of interest and proper background regions. Regardless of whether the background file is provided or not, HOMER will normalize the background to remove GC-bias and will also perform autonormalization. <br>

```
findMotifsGenome.pl [PROJECT-NAME-stageA-specific-TRE].bed \
hg38 \
[PROJECT-NAME-stageA-specific]_HOMER_output_folder \
-size given \
-bg [PROJECT-NAME-non-stageA-specific-TRE].bed \
-p 22
```

After the job is completed, find `homerResult.html` and `knownResults.html` in the output folder for result summary. You can find more informtion about how to inteprete the results on HOMER website, here are some key highlights: <br>
- **Known Motif Enrichment Results**: uses existing motif database (prior knowledge) to determine motif enrichment <br>
- **De Novo Motif Enrichment Results**: unbiased approach, does not take prior knowledge and finds repeating sequence segments in the data that are later compared with the known motif sequences <br>
- **% of Target Sequences with Motif**: % of motifs in stage of interest vs. custom background <br>
You can perform further filtering of HOMER results by p-value, q-value, fold enrichment, etc. to narrow motif list. <br>

### 6-3: Extract TREs containing specific motifs 
You can also extract TREs containing motif of "transcription factor X" using HOMER. Generage usage of the HOMER command is `annotatePeaks.pl <peak/BED file> <genome> -size # -m <motif>`. You will need to upload motif of interest from the output folder of HOMER motif enrichment analysis. <br> 

Below is the example for extracting stage A speicifc enhancers that contain motif of CDX2 in human samples. <br>
```
annotatePeaks.pl [PROJECT_NAME]_stageA_enhancer_coordinates.bed hg38 -size given -m CDX2.motif > HOMER_stageA_CDX2_enhancer.txt
```
This HOMER command generates a really detailed output file. We can use the python script to extract region coordinates for other downstream analyses. <br> 
```
python3 /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/extract_motif_coordinates_HOMER.py HOMER_stageA_CDX2_enhancer.txt > stageA_enhancer_containing_CDX2_motif.bed
```

### 6-4: Additional comments
In some cases were you have a specific hypothesis for a given TRE region (e.g. are there any motif sequences of transscript factor X present in a given TRE?), you can consider using function `FIMO` (Find Individual Motif Occurences) from [meme suite](https://meme-suite.org/meme/index.html). While `meme suite` has not been implemented in our lab's linux environment, you can easily carry out job using their website portal. <br>

You may need to do the following before execute FIMO job: <br>
- get DNA sequence of region of interest using this [tool](https://www.genome.ucsc.edu/cgi-bin/hgc?hgsid=1377211335_1JSdKvOBhA7HkslWuQydSAUvcXEN&g=getDna&i=NM_001328514.1&c=chr13&l=63278695&r=63795367&o=ncbiRefSeqCurated&table=ncbiRefSeqCurated) from UCSC Genome Browser. <br>
- download compatitle motif database for FIMO. In my case working with human samples, I downloaded `HOCOMOCOv11_full_HUMAN_mono_meme_format.meme` for motif scanning. <br>


## 7. TRE density analysis
>**Goal**: Count TREs Within X bp of each TSS/gene using GTF file as the genome reference. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/countTREsWithinXbp_v2.1.py` <br>
>**Appropriate type of compute node**: 24-core node <br>

We typically examine the TRE:TSS relationship by using search window within 100 kb upstream and downstream of TSS. 
An example command is:
```
python3 countTREsWithinXbp_v2.1.py [PROJECT_NAME]_stageA_enhancer_coordinates.bed \
/home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf \
100000 \
> [PROJECT_NAME]_stageA_enhancer_countTREsWithin100kb.txt
```


## 8. Super-enhancer analysis
>**Goal**: Define super-enhancers (or enhancer hotspots) with enhancers of interest. <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools_intragenicTRE/identifySuperEnhancers_intragenicTRE_v3.0.py` <br>
>**Appropriate type of compute node**: 24-core node <br>

### 8-1: Prepare files
1. Enhancers coordinates (bed format): see instructions in [HOMER section](#prepare-bed-files). <br>
2. Bigwig paths (txt format): line-delimited file with path to bigwig file for each sample of interest (per condition). Must include absolute **full** path to each bigwig file. Include only path to **plus** strand bigwig. Minus strand bigwig must have same prefix. <br>

### 8-2: Run super-enhancer (SE) analysis
This python script include the following major steps (*Note*: There are several differences compared with identifySuperEnhancers_v2.0.py): <br>
1. Input **enhancer** TREs of interest. The script doesn't take all TRE and remove promoter TREs for you. <br>
2. Stitch together TREs within defined distance (default: 12500 base) <br>
3. Get read counts within each TRE from bigwigs (count both strands for non-intragenic TREs; count the opposite strand and *2 for intragenic TREs) <br>
4. Normalizing read counts for differences in sequencing depth <br>
5. Sum TRE read counts within each stitched enhancer <br>
6. Rank enhancers and identify super enhancers <br>

Execute the the job by the following command: <br>
```
module load R/3.5.0
```
```
python3 /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools_intragenicTRE/identifySuperEnhancers_intragenicTRE_v3.0.py \
-p [PROJECT-NAME-SE-stageA-prefix] \
[PROJECT-NAME-stageA-specific-TRE].bed \
/home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf \
[PROJECT-NAME-path-to-stageA-bigwig].txt &> SE_output.log&
```

The completion of the job will generate `[PROJECT-NAME-SE-stageA-prefix]_rankedPlot.png` and `[PROJECT-NAME-SE-stageA-prefix]_stitched_enhancers_info.txt`. 

`[PROJECT-NAME-SE-stageA-prefix]_stitched_enhancers_info.txt` file contains the following columns: <br>
1. Parent: the coordinates of stitched enhancers 
2. Rank: the ranking based on the Total Signal of stitched enhancers
3. TotalSignal: the sum of the ChRO-seq signal from each of the individual enhancers within a stitched enhancer
4. Num TRE: the number of individual enhancers  stitched enhancer
5. Length: length of the parent coodinates 
6. Super: denotes whether the TotalSignal of a stitched enhancer is strong enough to be defined as a super-enhancer
7. Children: the coordinates of each individual enhancer of a given stitched enhancer


## 9. Assign genes to TREs/SEs
Assigning genes to TREs/SEs can be distance-based ([Method1](#9-2:-Method-1)) or correlation-based ([Method2](#9-3:-Method-2)). 

### 9-1: prepare bed files
Coordinates of TRE/SE/enhancer of interest in bed format: see instructions in [HOMER section](#prepare-bed-files).

### 9-2: Method 1 
>**Goal**: Identify genes of which the TSSs are closest to TREs/SEs of interest <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/findClosestGene2TRE_v2.0.py` <br>
>**Appropriate type of compute node**: 24-core node <br>

See an example command that assign genes to super-enhancers uniquely present in stage A. In the lab, we like to flag `-e (--express)` to run against a custome gene list instead of all genes across the entire genome. The gene list should be a line-delimited file containing a gene name per line. Depending on the purpose of this analysis could be genes that are transcribed in stage A (filtering by certain ChRO-seq signal threshold), are significantly up-transcribed in stage A compared to stage B (filtering based on DESeq2 analysis). Flag `-e` not only saves computing time but also generates more biological meaningful results.  
```
python3 /home/pr46_0001/cornell_tutorials/ChROseq_tutorial/findClosestGene2TRE_v2.0.py \
[PROJECT-NAME-stageA]_SE.bed \
/home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf \
-e [UP-genes-stageA-vs-stageB].txt \
> [PROJECT-NAME]_findClosestGene2TRE_stageA_SE_UPgene.txt
```

### 9-3: Method 2
>**Goal**: Performs Pearson correlation between genes and TREs within 500 kb (or specified distance) in distance following log2+1 transformation <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/GeneTRECorrelations2bedpe.py` <br>
>**Appropriate type of compute node**: 24-core node <br>

! Noted this correlation-based method need a large sample size to achieve statistical significance ! <br>

## 10. Define differentially transcribed (DT) genes between cell types/conditions
>**Goal**: Performs Pearson correlation between genes and TREs within 500 kb (or specified distance) in distance following log2+1 transformation <br>
>**Tool path**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/tools/gtf2bed_v2.1.py` <br>
>**R script template**: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ChROseq_R/DESeq2_ChROseq_Genes_tutorial.R` <br>
>**Appropriate type of compute node**: 24-core node <br>

### 10-1: Re-define gene body coordinates
RNA polymerase pausing around promoter regions has been observed (see [Core at al., Science, 2008](https://www.science.org/doi/10.1126/science.1150843)). To eliminate the bias generated from this RNA pol pausing event, we exclude ChRO-seq signal at the promoter regions when extracting counts from gene bodies. Using a python script, we re-define gene bodies coordinates by removing the first 500 bp (or custom distance) downstream of all TSSs in gtf annotations. An example command will be: <br>
```
python3 gtf2bed_v2.1.py -s 500 /home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.annotation.gtf > [PROJECT_NAME]_gencode_annotation_start500.bed
```
You can also flag the python script to remove ChRO-seq signal around the transcription terminal sties (TTSs). However, based on my ChRO-seq datasets, I find the RNA pol pausing at TTSs is not profound. The output bed file is required for the R script described in the next session. <br>

### 10-2: Define DT genes between cell types/conditions
The R script template of this job can be found at: `/home/pr46_0001/cornell_tutorials/ChROseq_tutorial/ChROseq_R/DESeq2_ChROseq_Genes_tutorial.R`. <br>

See detailed instructions [here](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=266#c) for using Rstudio at Cornell BioHPC. Make sure to load R 3.5.0 for this job. <br>
```
/programs/rstudio_server/rstudio_start 3.5.0
```

You can modify the R template for your project needs. The R script include the 3 major steps (see below). <br>
1. Extract ChRO-seq counts from gene bodies using package `bigwig` <br>
2. Exploratory data analysis including heirachical clustering analysis and data dimentionality reduction by PCA 
3. Input count matrix for differential expression analysis using package `DESeq2` <br>

### 10-3: Additional comments
In addition to RNA pol pausing at TSS/TTS regions, the active TREs present in gene bodies may create biases when extracting counts from gene bodies. One way to eliminate the bias from the activity of intragenic TREs is to exclude counts from these regions. I have done this for [Hung et al., bioRxiv (multi-omics project)](https://www.biorxiv.org/content/10.1101/2022.07.12.499825v1.full). Below are some technical details about how to acheive it: <br>

1. Use `gtf2bed.py` to re-define the gene coordinates (ex: set gene start position as 500 bp downstream of TSS) <br>
2. Remove new gene coordinates with negative length. A simple `awk` command can tackle this: <br>
```
awk '$5 > 0 {print $0}' genecoordinates_start500.bed > genecoordinates_start500_positive.bed
```
3. Subtract TRE regions from gene bodies. A simple `bedtool` command can tackle this: <br>
```
bedtools subtract -a genecoordinates_start500.positive.bed -b TREs_coordinates.bed > genecoordinates_start500_positive_noTRE.bed
```
Say if you want to compare stage A relevant to stage B, you would revmoe TERs that are present in one or both stages. 


## 11. Other ChRO-seq related tools

The [Danko lab](https://www.dankolab.org/Software.html) is very active in expanding and maintaining analysis tools for ChRO-seq downstream analysis, so you can check it out if you like! The one that 

### Tunits

Classify predicted transcription units: Haven't polished python script for publication purpose or 

