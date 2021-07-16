# GeneAnnotations (GGA)
This script will define active transcripts and annotate dominant TSS and TES locations using RNA-seq and PRO-seq/Start-seq data. 

## General workflow
To call active transcripts and dominant TSS and TES using expression data, this script will:

1) Use RNA-seq data to quantify transcript expression, using Kallisto to quantify reads over individual transcripts.

2) If  5' end PRO-seq or Start-seq data is available, use TSScall to generate a list of active TSSs and then determine the dominant TSS per gene (non-dominant TSS information is retained). 
- If you only have RNA-seq data, your transcript TPM values will be used to call dominant TSSs.

3) For the dominant TSS and associated transcripts, the script uses the average transcript TPM values to select the dominant TES.

The end result is the designation of a single, dominant TSS and TES window for each expressed gene, giving a 1:1 correspondence between genes active in RNA-seq and genomic TSS and TES locations to investigate using PRO-seq, ChIP-seq, etc. 

Please note: 
- The pipeline produces files in the working directory, so we recommend running it from a folder on scratch. The pipeline generates an extensive number of temporary files during runtime. Please ignore these files. 
- When the script is complete, there will be 5 final files in your working directory. Please refer to "Pipeline outputs and critical information for downstream analysis" for details on each of these files. 
- The pipeline runs a series of sbatch scripts, most of which won't start until the previous one has successfully finished. This means that you will see a lot of jobs queued with the "(Dependency)" listed. In the **logs** folder, you can find detailed log files for each step in the pipeline. The number at the end of the log filename corresponds to the slurm job id. 
 
## Running the pipeline

To run the gene annotations pipeline you need a few things

**1) A manifest text file of RNA-seq data:**
- 3 tab-delimited columns, 1 row per RNA-seq sample
- 1st column: full path to R1 fastq
- 2nd column: full path to R2 fastq
*Note: multiple fastqs per sample can be given as comma-separated list (eg. if the same sample was sequenced on multiple runs), just ensure that R1 and R2 fastq files are listed in the same order.
*Note: For single-end RNA-seq, simply input a 2-column manifest (fastq and samplename columns). 
- 3rd column: unique sample name

**2) 5' mapped PRO-seq fwd and reverse bedgraphs (no normalization, all relevant files concatenated together, single bp resolution)**
IMPORTANT: Please refer to PRO-seq/TSScall FAQ #2. For 5'end mapping, it is required that PRO-seq libraries were sequenced as paired-end.

Note: We have observed that TSScall can time-out if the input F and R bedgraphs are extremely large. To resolve this issue, we recommend that you remove positions from the input bedGraphs that have very few reads. The number you will use for filtering bedgraphs will depend on the number of samples merged and the sequencing depth. As a starting benchmark, we recommend removing positions with fewer than 4 reads (Please see example below). If TSScall still times-out, then we recommend removing positions with fewer than 7 reads (The default read threshold per position in TSScall is 8). 

```
awk '{if ($4 >= 4) print $0}' orig_F.bedGraph > filtered_F.bedGraph
awk '{if ($4 >= 4) print $0}' orig_R.bedGraph > filtered_R.bedGraph
```

**3) Pipeline input arguments (fastq trim positions, reference genome build, and output prefix)**

The full list of options and usage is below:
```
usage: get_gene_annotations.sh [-h | --help]
                  [-m | --manifest MANIFEST_FILE]
                  [-R | --refgenome REF_GENOME_BUILD]
                  [-o | --output OUTPUT_PREFIX]
                  [-x | --mintrim MIN_FASTQ_TRIM_POSITION]
                  [-X | --maxtrim MAX_FASTQ_TRIM_POSITION]
                  [-t | --tpm RNA-SEQ_TPM_THRESHOLD]
                  [-n | --num NUMBER_OF_SAMPLES_NEEDING_TO_PASS_TPM_THRESHOLD]
                  [-f | --proseqFwd 5PRIME_PROSEQ_FWD_BEDGRAPH]
                  [-r | --proseqRev 5PRIME_PROSEQ_REV_BEDGRAPH]
                  [--skipfastqfilter]
                  [--noProSeq]
                  [--stopAtGTFfilter]
                  [--single]

required arguments:
  -m, --manifest        tab-delimited text file with R1 fastq, R2 fastq, and samplename.
                        Note: RNA-seq Fastq files can be .gz or .bz2 zipped. can give multiple
                        fastq files per R1 and R2 as comma-separated list (but they must
                        be in same order for both R1 and R2 lists)
  -R, --refgenome       reference genome build (mm10|hg19|hg38)
  -o, --output          output prefix for created files
  -x, --mintrim         leftmost trim position for fastq filter (usually 1) -- not required when using --single
  -X, --maxtrim         rightmost trim position for fastq filter (usually read length) -- not required when using --single
  -f, --proseqFwd       5' mapped PRO-seq forward strand bedgraph (filtered, e.g. > 4 reads)
  -r, --proseqRev       5' mapped PRO-seq reverse strand bedgraph (filtered, e.g. > 4 reads)


optional arguments:
  -h, --help            show this help message and exit
  -t, --tpm             TPM threshold for transcript filtering (default = 0)
  -M, --tsscallMin      TSScall minimum read count (default = 8)
  -w, --tsscallWindow   TSS seachh window in bp (default = 1000)
  -T, --tsscallTime     Time (in minutes) for TSScall sbatch job, default = 1500
  -P, --tsscallPartition Partition for TSScall sbatch job, default = medium
  -s, --stranded        Stranded option for RNA-seq mapping by Kallisto. Choose one of 
                        fr-stranded, rf-stranded, or unstranded (default = rf-stranded)
  --skipfastqfilter     skip filtering step for fastq, note that fastqs files must be
                        unzipped to use this option
  --noProSeq            If you only have RNA-seq, use this option to skip the PRO-seq
                        steps
  --stopAtGTFfilter     Stop pipeline after gtf filtering, before TSScall
  --single              Single End RNA-seq fastqs. Doesn't trim/filter, just proceeds directly to kallisto. Manifest should just have two columns (fastq and samplename). 

```

### Pipeline example run

Run Pipeline: 
```
{path}/GGA.sh --manifest metadata.txt --refgenome mm10 --output Name -x 1 -X 75 --tpm 0 --num 1 --proseqFwd {path}/PROseq_5readmap_F.bedGraph --proseqRev {path}/PROseq_5readmap_R.bedGraph
```
Note, for single-end data make a manifest with just two columns (fastq and samplename) and specify the --single option for the *get_gene_annotations.sh* script. Reminder: If you are using single-end RNA-seq data, you will need to provide trimmed, filtered and unzipped fastqs. 

## Dependencies

The pipeline is set up to run in the HMS O2 servers, using the pre-compiled modules and Slurm commands of that environment. To run these scripts outside of this environment please ensure you have the appropriate dependencies available and ensure the Slurm sbatch commands are properly configured. Also make sure that the Stringr package is installed in R. 

The pipeline also requires the transcript GTF and kallisto indexes to already be made and in the correct directory. See the **How_to_make_cDNA_fasta_files.md** markdown file for more info on how these are made. The paths to the scripts and indexes should be specified in *get_gene_annotations.sh*

## Pipeline outputs and critical information for downstream analysis 

The output files all make use of the *<output_prefix>* parameter passed to the pipeline script. The tables below describes the output files

1. _**<output_prefix>_Dominant.TSS.TES.calls_WORKING_FILE.txt**_ : Dominant TSS and TES information for each gene supported by RNA-seq and PRO-seq (if provided). **THIS IS THE MAIN WORKING FILE!!** For both minus and plus strand annotations, the dominaint TSS column refers to the TSS coordinate and the dominant TES column refers to TES coordinate. In all other output formats (bed, and txt for makeheatmap), the coordinates are listed as min and max. This means that for plus strand annotations the min column refers to the dominant TSS, whereas for - strand annotations the min coordinate refers to the dominant TES. **To prevent confusion, we recommend that ANY additional extraction of dominant TSS and TES coordinates should ONLY come from this file.**
2. _**<output_prefix>_Dominant.TSS.TES.calls.bed**_ : BED file of <output_prefix>_Dominant.TSS.TES.calls_WORKING_FILE.txt. This file can be used to upload the dominant annotation to Genome Browser. Column 2 indicates the min coordinate and Column 3 indicates the max coordinate. 
3. _**<output_prefix>_Dominant.TSS.TES.calls_formakeheatmap.txt**_ : TSS to TES file for makeheatmap (Column 2 anchor = TSS, Column 3 = Min coordinate, Column 4 = Max coordinate).
4. _**<output_prefix>_Dominant.Affiliated.ActiveTranscripts.gtf**_ : GTF for all transcripts stemming from the dominanat TSS to TES (can be multiple transcripts per gene). This GTF can be used for featurecounts
5. _**<output_prefix>_Annotated_Dominant_and_Nondominant_obsTSS_fordREG.txt**_ : List of all active and annotated obsTSSs identified by TSScall. Each obsTSS was annotated for (A) If the obsTSS was selected as the dominant TSS per Gene (Dominant_obsTSS Column = TRUE if Dominant TSS) and (B) If the obsTSS is associated with stable and active transcripts defined by RNA-seq (Stable Column = TRUE if associated with a stable and active transcript). This file can be used to overlap active and annotated obsTSSs with features such as dREG or ATAC-seq peaks. Please note that both stable (associated transcript(s) supported by RNA-seq) and unstable (associated transcript(s) not supported by RNA-seq) gene models are included in this list. The final number of dominant obsTSSs pulled from this file will not match the number from <output_prefix>_Dominant.TSS.TES.calls_WORKING_FILE.txt

### General FAQs

#### What is the appropriate GTF to use for featurecounts?
Our script will output a GTF that contains active transcript models that stem from the dominant TSS to TES clusters. 

#### How to quantify PRO-seq read counts across genes?
The first step for PRO-seq analysis is to count sense-strand reads from TSS to TES. This quantification is best done using makeheatmap (see its readme here: __https://github.com/AdelmanLab/NIH_scripts/tree/main/make_heatmap__) to count reads from your raw 3' PRO-seq bedGraphs. For this quantification:
1. Use the _**<output_prefix>_Dominant.TSS.TES.calls_formakeheatmap.txt**_ file generated by GGA.
1. Quantify reads in each sample using makeheatmap: run in variable bin size mode, specifying just 1 bin, and counting reads on the sense strand only.
1. After running makeheatmap on each sample, paste the counts from each file together into a read count matrix. This can then be input into DESeq2 for differentially expressed gene analysis.

#### What are all the Columns in the TSScall matrix file?: _**<output_prefix>_Merged.FINAL.Dominant.TSS.TES.calls.txt**_
1. *GeneID:* The ENSEMBL Gene ID
2. *GeneName:* The Gene Symbol (or name) for each ID
3. *GeneType:* The annotated gene biotype
4. *TSS_ID:* The TSScall ID for the dominant TSS
5. *Chr:* Chromosome
6. *Strand:* Plus or minus strand
7. *TSScallOnlyUnstableTx:* Flag for whether the transcripts stemming from the dominanat TSS to TES are not associated with RNA-seq counts (above the specificed threshold) and are presumably unstable (1 if it is unstable). In the script outputs all entries should have a zero in this column. 
8. *DominantTSS:* The position of the Dominant TSS 
9. *DominantTES:* The position of the Dominant TES 
10. *TSScallReads:* The number of 5' PRO-seq reads at the TSS position. Note* this number refers specifically to the TSS nucleotide. This is NOT the number you will get when summing reads around the TSS window (for example, from +/- 100nt).
11. *Distance.kallisto.to.TSScall.TSS:* Distance from the TSScall TSS to the annotated TSS identified by Kallisto using the RNA-seq and annotations
12. *Divergent:* Flag for whether there is a divergent TSS (TRUE/FALSE)
13. *Divergent.partner:* TSS call ID for the divergent TSS
14. *Divergent.distance:* Distance to the divergent TSS
15. *Convergent:* Flag for whether there is a convergent TSS (TRUE/FALSE)
16. *Convergent.partner:* TSS call ID for the convergent TSS
17. *Convergent.distance:* Distance to the convergent TSS
18. *MultipleGeneID:* Flag for whether the TSS is associated with multiple Gene IDs (TRUE/FALSE)
19. *AlternativeGeneID:* Gene IDs also associated with dominant TSS (if any)
20. *TESdominantTxITESclusterTPM:* RNA-seq TPM score for the dominant TES cluster
21. *TESallTxIDs:* All transcript models starting from the dominant TSS cluster and ending at the dominant TES cluster. 

## What genome Annotations are used?
By default the pipeline starts with the ENSEMBL genecode basic annotations. 

If you need to generate your own GTFs and cDNA fasta files, see the **How_to_make_cDNA_fasta_files.md** markdown file for example code.


## What FASTQ files are appropriate for GGA?
By default the FASTQ files are trimmed and filtered before running kallisto. For more on this check out the **fastq_process.sh** script. If you are using single-end RNA-seq data, you will need to provide trimmed, filtered and unzipped fastqs. 

## PRO-seq/TSScall FAQ: 

**Q1. What is TSScall?** 

TSScall was originally developed by Andy Lavender to idenify transcription start sites (TSS) from Start-seq data. TSScall uses a user-defined search window to identify active annotated TSSs (obsTSS) from reference annotations.

Please refer to the TSScall github for more specifics on TSScall: https://github.com/lavenderca/TSScall

**Q2. If TSScall was designed for Start-seq, how can PRO-seq datasets be applied to call annotated active TSSs with TSScall?** 

Normally PRO-seq bedgraphs report the 3'end position of a read to specify the position of RNAPII. However, since the short promoter proximal RNAs are largely unaffected by base hydrolysis, PRO-seq data where the 5'end of the read is reported can be used to accurately identify obsTSSs when Start-seq datasets are not available. IMPORTANT: For 5'end mapping, it is required that PRO-seq libraries were sequenced as paired-end. Additionally, TSScall using PRO-seq 5'end read information can be used to identify convergent and divergent obsTSSs. 
 
**Q3. What is the recommended read count threshold?** 

We recommend setting set_read_threshold as 8, however this number can be lowered based on available sequencing depth. 

**Q4. Can TSScall identify divergent and convergent TSSs?** 

Yes! In order to identify the maximum number of uaTSSs and convergent TSSs, we recommend running TSScall with a low set_read_threshold (default threshold = 8 reads), and then stringently filter obsTSSs after TSScall has been run. 

**Q5. Can I trust nuTSS annotations called from PRO-seq data?**  

nuTSS annotations are ignored when TSScall is run on PRO-seq data because PRO-seq 5' ends within gene bodies can inadvertently be identified by TSScall as nuTSSs. If a nuTSS is a true uaTSS or convergent TSS it will be accurately paired with an obsTSS, and since these fall on the antisense strand, are not subject to the same problems as sense-strand calls. Information on the position of these divergent and convergent TSSs is included in the TSScall output.

## Kallisto FAQ

Kallisto uses pseudoalignmnet for a probabilistic quantification of transcript abundance. This approach makes it both fast and accurate. For more on Kallisto, see the about section of the kallisto webpage: https://pachterlab.github.io/kallisto/about

The pipeline uses kallisto to quantify transcript abundance in each of the submitted RNA-seq samples. Then for TES calling, the average transcript TPM across all the samples is used. Starting from the dominant TSS associated transcripts, the annotated transcript end sites within 100bp of each other are clustered. The TPM values for clustered transcripts are added together, and the cluster with the highest TPM is selected as the dominant cluster. For clusters with multiple transcripts, the transcript with the highest TPM is selected as the dominant transcript and the TES for that transcript is reported as the dominant TES.

For RNA-seq derived TSS calling, the average transcript TPM across all the samples is used. The annotated transcript start sites within 100bp of each other are clustered. The TPM values for clustered transcripts are added together, and the cluster with the highest TPM is selected as the dominant cluster. For clusters with multiple transcripts, the transcript with the highest TPM is selected as the dominant TSS. The transcripts in the dominant TSS cluster are then passed on to TES calling as above (^).






