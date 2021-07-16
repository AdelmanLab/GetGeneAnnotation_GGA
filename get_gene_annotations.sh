#!/bin/bash
### this is not updated (obviously), once script is finished will update this. this runs if error with inputs or help option specified
usage()
{
    echo "usage: get_gene_annotations.sh [-h | --help]
                  [-m | --manifest MANIFEST_FILE]
                  [-R | --refgenome REF_GENOME_BUILD]
                  [-o | --output OUTPUT_PREFIX]
                  [-x | --mintrim MIN_FASTQ_TRIM_POSITION]
                  [-X | --maxtrim MAX_FASTQ_TRIM_POSITION]
                  [-t | --tpm RNA-SEQ_TPM_THRESHOLD]
                  [-n | --num NUMBER_OF_SAMPLES_NEEDING_TO_PASS_TPM_THRESHOLD]
                  [-f | --proseqFwd 5PRIME_PROSEQ_FWD_BEDGRAPH]
                  [-r | --proseqRev 5PRIME_PROSEQ_REV_BEDGRAPH]
                  [-T | --tsscallTime TSSCALL_SBATCH_TIME] 
                  [-P | --tsscallPartition TSSCALL_SBATCH_PARTITION] 
                  [--skipfastqfilter]
                  [--noProSeq]
                  [--stopAtGTFfilter]

required arguments:
  -m, --manifest        tab-delimited text file with R1 fastq, R2 fastq, and samplename.
                        Note, Fastq files can be .gz or .bz2 zipped. can give multiple
                        fastq files per R1 and R2 as comma-separated list (but they must
                        be in same order for both R1 and R2 lists)
  -R, --refgenome       reference genome build (mm10|hg19|hg38)
  -o, --output          output prefix for created files
  -x, --mintrim         leftmost trim position for fastq filter (usually 1)
  -X, --maxtrim         rightmost trim position for fastq filter (usually read length)
  -f, --proseqFwd       5' mapped PRO-seq forward stand bedgraph
  -r, --proseqRev       5' mapped PRO-seq reverse stand bedgraph


optional arguments:
  -h, --help            show this help message and exit
  -t, --tpm             TPM threshold for transcript filtering (default = 0)
  -n, --num             The number of samples needing to pass TPM threshold (default = 1)
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
    "
}

#set scripts and annotations paths
#note, these the paths for these folders on your system will have to be specified for the script to run
   scriptsPath=''
   annotationPath=''
   rlibpath=''
   chrsizePath=${scriptsPath}/chromSizes



################### Read in options


#some of these are setup as defaults, which can be overriden, while others need to be entered

manifest=
gtf=
processfastq=1
stranded=rf-stranded
single=0
expname=
ref=
kalIndex=
tpm=0
numSamples=1
mintrim=
maxtrim=
proseqFwd=
proseqRev=
tsscallMin=8
tsscallWindow=1000
useProseq=1
stopAtGTFfilter=0
tsscallTime=1500
tsscallPartition=medium
gtf=
kalIndex=
txInfo=
chrsizes=
admin=0

while [ "$1" != "" ]; do
    case $1 in
         -m | --manifest )     shift
                                manifest=`readlink -e $1`
                                ;;
#        -g | --gtf )           shift 
#                                gtf=$1
#                                ;;
        -R | --refgenome )      shift
                                ref=$1
                                ;;
        -t | --tpm )           shift
                                tpm=$1
                                ;;
        -n | --num )           shift
                                numSamples=$1
                                ;;
        -s | --stranded )      shift
                                stranded=$1
                                ;;
        -o | --output )        shift
                                expname=$1
                                ;;
        -x | --mintrim )       shift
                                mintrim=$1
                                ;;
        -X | --maxtrim )       shift
                                maxtrim=$1
                                ;;
        -f | --proseqFwd )     shift
                                proseqFwd=$1
                                ;;
        -r | --proseqRev )     shift
                                proseqRev=$1
                                ;;
        -M | --tsscallMin )    shift
                                tsscallMin=$1
                                ;;
        -w | --tsscallWindow ) shift
                                tsscallWindow=$1
                                ;;
        -T | --tsscallTime ) shift
                                tsscallTime=$1
                                ;;
        -P | --tsscallPartition ) shift
                                tsscallPartition=$1
                                ;;
        --gtf )                shift
                                gtf=$1
                                ;;
        --kalIndex )           shift
                                kalIndex=$1
                                ;;
        --txInfo )             shift
                                txInfo=$1
                                ;;
        --chrsizes )           shift
                                chrsizes=$1
                                ;;
        --admin )                shift
                                admin=$1
                                ;;
        --single )             shift
                                single=1
                                ;;
        --stopAtGTFfilter )    shift
                                stopAtGTFfilter=1
                                ;;
        --skipfastqfilter )         shift
                                processfastq=0
                                ;;
        --noProSeq )         shift
                                useProseq=0
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done




echo "reading in options"

#if no gtf given, then set gtf, etc, by ref genome
if [ -z "$gtf" ];then
   gtf=${annotationPath}/${ref}/${ref}.basic.gtf
   kalIndex=${annotationPath}/${ref}/${ref}.basic.cDNA.idx
   txInfo=${annotationPath}/${ref}/${ref}.basic.transcripts.to.genes.info.txt
   chrsizes=${chrsizePath}/${ref}.chrom.sizes
   cleanChrsizes=${chrsizePath}/Clean/${ref}.chrom.sizes
fi


#kallisto strandedness options

if [ "$stranded" == "rf-stranded" ]; then
   kalStrand='--rf-stranded'
elif [ "$stranded" == "fr-stranded" ]; then
   kalStrand='--fr-stranded'
elif [ "$stranded" == "unstranded" ]; then
   kalStrand=''
fi




echo "ref genome: $ref"
echo "gtf: $gtf"
echo "kallisto index: $kalIndex"
echo "transcript info file; $txInfo"

echo "manifest file: $manifest"
echo "min trim: $mintrim"
echo "max trim: $maxtrim"

if [ "$single" -eq 1 ]; then
   echo "single-end RNA-seq option selected"
fi
if [ "$single" -eq 0 ]; then
   echo "paired-end RNA-seq option selected"
fi

mkdir -p logs


echo "loading modules"
module purge
module load gcc/6.2.0
module load python/2.7.12
module load cutadapt/1.14
module load kallisto/0.45.1
module load R/4.0.1
module load python/2.7.12


#clean up manifest file


manifest=`realpath $manifest`

endline=`tail -c 1 "$manifest"`
if [ -z "$(tail -c 1 "$manifest")" ]; then
        echo "Newline already at end of file"
    else
        echo "Adding newline at end of file"
        echo -e "\n" >> $manifest
    fi

#make expname folder
mkdir -p ${expname}


#added in an option for single-end RNA-seq data 
# here is the paired-end processing, this if statement runs through the end of the kallisto quant steps
if [ "$single" -eq 0 ]; then


################### Process fastqs (Paired end)
## This step can be skipped by specifying --skipfilter 
if [ "$processfastq" -eq 1 ]; then
   mkdir -p ${expname}/fastq
   cd ${expname}/fastq
   
   echo "filtering and trimming fastqs"

          {
                while read -r f1 f2 f3
                do
                     if [[ -n "$f3" ]]; then #only submit jobs with sample name
                        fileR1=$f1
                        fileR2=$f2
                        samplename=$f3
                        echo $samplename
                        #assemble a list of the jobids for the sbatch submissions (kallisto_jobids)
                        fastq_jobids=""
                        kjob_id=($(sbatch -p short --mem=30G -t 360 --parsable --output=../../logs/${samplename}_fastq_process_%j.log ${scriptsPath}/fastq_process_and_kallisto_map.sh -s $scriptsPath --read1 $fileR1 --read2 $fileR2 -m $mintrim -M $maxtrim --outprefix $samplename --stranded ${stranded} --kalIndex ${kalIndex}))
                        if [ -z $kallisto_jobids ];then
                                kallisto_jobids=`echo "${kjob_id}"`
                                echo "first jobid"
                        else
                                kallisto_jobids=`echo "${kallisto_jobids},${kjob_id}"`
                                echo "added job id to list"
                        fi
                   fi
                done
        } < "$manifest"

   cd ../..
fi

fi


################### run kallisto on each sample
#note this is a bit ugly... i have an if/else statment for whether the fastq files were trimmed/filtered to begin with or not. this changes how i find the fastq files and whether the sbatch is dependent on the previous step or not.
if [ "$single" -eq 0 ]; then
   if [ "$processfastq" -eq 0 ]; then

   cd ${expname}

echo -e "\n\n setting up paired-end kallisto sbatches"
{
while read -r f1 f2 f3
                do
                     if [[ -n "$f3" ]]; then #only submit jobs with sample name
                        
                        samplename=$f3                
                        untrimmedR1=$f1
                        untrimmedR2=$f2
                        kjob_id=($(sbatch --parsable --output=../logs/${samplename}_kallisto_%j.log  -t 120 -p short --mem=15G -c 5 --wrap "kallisto quant ${kalStrand} -i ${kalIndex} -o $samplename -t 5 -b 30 ${untrimmedR1} ${untrimmedR2}"))
                        if [ -z $kallisto_jobids ];then
                                kallisto_jobids=`echo "${kjob_id}"`
                                echo "first jobid ${kjob_id}"
                        else
                                kallisto_jobids=`echo "${kallisto_jobids},${kjob_id}"`
                                echo "added ${kjob_id} job id to list"
                        fi
                        
                     fi
                done
        } < "$manifest"

   cd ..

   fi
fi

################### Kallisto on single end data

if [ "$single" -eq 1 ]; then

cd ${expname}

echo $manifest

echo -e "\n\n setting up single-end kallisto sbatches"
{
while read -r f1 f2
                do
                     if [[ -n "$f2" ]]; then #only submit jobs with sample name
                        fastq=$f1
                        samplename=$f2
                        ls $fastq
                        echo $samplename
                        
                        fragmentlength=200
                        fragmentstdev=30
                        
                        kjob_id=($(sbatch --parsable --output=../logs/${samplename}_kallisto_%j.log  -t 120 -p short --mem=15G -c 5 --wrap "kallisto quant ${kalStrand} -l ${fragmentlength} -s ${fragmentstdev} -i ${kalIndex} -o $samplename --single -t 5 -b 30 ${fastq}"))
                        if [ -z $kallisto_jobids ];then
                                kallisto_jobids=`echo "${kjob_id}"`
                                echo "first jobid ${kjob_id}"
                        else
                                kallisto_jobids=`echo "${kallisto_jobids},${kjob_id}"`
                                echo "added ${kjob_id} job id to list"
                        fi
                        
                     
                   fi
                done
        } < "$manifest"

cd ..

fi



################### filter GTF with kallisto counts
echo "kallisto job ids: ${kallisto_jobids}"

#run Rscript to filter GTF and get Avg TPM values per transcript
### note, now removing weird chromosomes at this step
echo -e "\n\n setting up GTF filter script sbatch"

echo -e '#!/bin/bash\n#SBATCH -p short\n#SBATCH -t 300\n#SBATCH --mem=20G\n#SBATCH \n#SBATCH -o r_scripts.out\n\nRscript '${scriptsPath}'/gtf_filter_from_tpms.R '${tpm}' '${numSamples}' '${gtf}' '${expname}' '${txInfo}' '${cleanChrsizes}'' > gtf_filter_script.sh

gtf_filter_job_id=($(sbatch --output=logs/${expname}_gtf_filter_%j.log --parsable --dependency=afterok:$kallisto_jobids gtf_filter_script.sh))



#stop here if stopAtGTFfilter option given
if [ "$stopAtGTFfilter" -eq 1 ]; then
   exit 0
fi


############################## call dominant TSS and TES solely from RNA-seq, this step also filters to only keep annotations on the cleaned up chromosome list (redundant but keeping it in here as well)

echo "call TSS and TES dominant from RNA-seq kallisto counts"

echo -e '#!/bin/bash\n#SBATCH -p short\n#SBATCH -t 150\n#SBATCH --mem=20G\n#SBATCH \n#SBATCH -o r_scripts.out\n\nRscript '${scriptsPath}'/kallisto_tss_tes_call.r '${expname}'_avg_tpms_all_tx.txt '${gtf}' '${expname}' '${cleanChrsizes}'' > kallisto_tss_tes_rscript.sh

kallist_tss_tes_call_job_id=($(sbatch --output=logs/${expname}_kallist_tss_tes_call_rscript_%j.log --parsable --dependency=afterok:$gtf_filter_job_id kallisto_tss_tes_rscript.sh))





##################### TSS Call and clean up
if [ "$useProseq" -eq 1 ]; then

echo "setting up TSScall jobs"

tsscall_jobid=($(sbatch -p ${tsscallPartition} -t ${tsscallTime} --mem=80G --output=logs/${expname}_TSScall_%j.log --parsable --dependency=afterok:$kallist_tss_tes_call_job_id --wrap "python ${scriptsPath}/TSScall-master/TSScall.py --annotation_file ${gtf} --detail_file ${expname}_TSScall_detail_file --set_read_threshold ${tsscallMin} --annotation_join_distance 500 --annotation_search_window ${tsscallWindow} ${proseqFwd} ${proseqRev} ${chrsizes} ${expname}_TSScall_output.bed"))

tssclassify_jobid=($(sbatch -p short -t 300 --mem=50G --output=logs/${expname}_tssclassify_%j.log --parsable --dependency=afterok:$tsscall_jobid --wrap "perl ${scriptsPath}/TSScall-master/TSSclassify.pl ${expname}_TSScall_detail_file ${gtf} > ${expname}_TSScall_detail_file_TSSclassify"))

## note this step removes the weird chromosomes from the TSS call outputs
echo -e '#!/bin/bash\n#SBATCH -p short\n#SBATCH -t 300\n#SBATCH --mem=20G\n#SBATCH \n#SBATCH -o r_scripts.out\n\nRscript '${scriptsPath}'/TSScall_cleanup_script.r '${expname}'_TSScall_detail_file_TSSclassify '${expname}' '${rlibpath}' '${cleanChrsizes}' '${tpm}'' > tsscall_cleanup_rscript.sh

tss_clean_dominant_job_id=($(sbatch --output=logs/${expname}_tsscall_cleanup_dominant_rscript_%j.log --parsable --dependency=afterok:$tssclassify_jobid tsscall_cleanup_rscript.sh))




#### Use kallisto counts to call dominant TESs on TSS call

echo -e '#!/bin/bash\n#SBATCH -p short\n#SBATCH -t 150\n#SBATCH --mem=20G\n#SBATCH \n#SBATCH -o r_scripts.out\n\nRscript '${scriptsPath}'/tsscall_kallisto_tes.r '${expname}'_avg_tpms_all_tx.txt '${expname}'_FINAL_Dominant_obsTSS_perGeneID_List.txt '${expname}' '${gtf}'' > tsscall_TES_rscript.sh

tsscall_TES_call_job_id=($(sbatch --output=logs/${expname}_tsscall_TES_call_rscript_%j.log --parsable --dependency=afterok:$tss_clean_dominant_job_id tsscall_TES_rscript.sh))



fi





############################### Merge TSScall and Kallisto derived annotations

if [ "$useProseq" -eq 1 ]; then
  tsscall_TES_call_kallist_tss_tes_call_job_ids=`echo "${tsscall_TES_call_job_id},${kallist_tss_tes_call_job_id}"`
  
  echo -e '#!/bin/bash\n#SBATCH -p short\n#SBATCH -t 150\n#SBATCH --mem=20G\n#SBATCH \n#SBATCH -o r_scripts.out\n\nRscript '${scriptsPath}'/TSScall_kallisto_meged_annotations.r '${expname}'_Kallisto_called_dominant_TSS_to_TES_clusters.txt '${expname}'_FINAL_Dominant_obsTSS_and_TES_perGeneID_List.txt '${expname}'_avg_tpms_all_tx.txt '${expname}' '${tpm}'' > kallisto_tsscall_merge_rscript.sh

  merge_job_id=($(sbatch -p short -t 30 --mem=10G --output=logs/${expname}_TSScall_kallist_merge_%j.log --parsable --dependency=afterok:$tsscall_TES_call_kallist_tss_tes_call_job_ids kallisto_tsscall_merge_rscript.sh ${expname}))


#Rscript ${scriptsPath}/TSScall_kallisto_meged_annotations.r ${expname}_Kallisto_called_dominant_TSS_to_TES_clusters.txt ${expname}_FINAL_Dominant_obsTSS_and_TES_perGeneID_List.txt ${expname}/${expname}_avg_tpms_all_tx.txt ${expname}

# 1) Kallisto TSS/TES calls (use the filtered one here)
# 2) TSScall final output
# 3) Kallisto Tx counts all
# 3) output prefix
  

fi




##### get GTF for TSS to TES transcripts
if [ "$useProseq" -eq 1 ]; then
  
  echo -e '#!/bin/bash\n#SBATCH -p short\n#SBATCH -t 45\n#SBATCH --mem=20G\n#SBATCH \n#SBATCH -o r_scripts.out\n\nRscript '${scriptsPath}'/get_dominant_TSS_to_TES_stable_transcripts_gtf.r '${expname}'_Dominant.TSS.TES.calls_WORKING_FILE.txt '$gtf' '${expname}'_Dominant.Affiliated.ActiveTranscripts_forfeaturecounts.gtf' > make_tss_to_tes_transcripts_gtf_rscript.sh

  tss_tes_gtf_job_id=($(sbatch -p short -t 30 --mem=10G --output=logs/${expname}_make_tss_to_tes_transcripts_gtf_%j.log --parsable --dependency=afterok:$merge_job_id make_tss_to_tes_transcripts_gtf_rscript.sh ${expname}))

fi



############################## Clean up files at the end

if [ "$useProseq" -eq 1 ]; then
 
  final_file_cleanup_job_id=($(sbatch -p short -t 20 --mem=1G --output=logs/${expname}_final_file_cleanup_%j.log --parsable --dependency=afterok:$tss_tes_gtf_job_id ${scriptsPath}/clean_up_files.sh ${expname} ${admin}))

else
  final_file_cleanup_job_id=($(sbatch -p short -t 20 --mem=1G --output=logs/${expname}_final_file_cleanup_%j.log --parsable --dependency=afterok:$kallist_tss_tes_call_job_id ${scriptsPath}/clean_up_files.sh ${expname} 1))

fi


