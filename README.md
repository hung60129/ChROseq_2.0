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
>**Goal**: Define a universal set of transcriptional regulatory elements (TREs) using the latest [dREG tool](https://github.com/Danko-Lab/dREG/) from the Danko lab.  <br>
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
### Greate a genome track on UCSC genome browser
### Other tools


## 5. TRE type classification


## 6. Define differentially transcribed TREs between cell types/conditions


## 7. Super-enhancer analysis


## 8. Motif enrichment analysis


## 9. Define differentially transcribed genes between cell types/conditions


## 10. Other ChRO-seq related tools


