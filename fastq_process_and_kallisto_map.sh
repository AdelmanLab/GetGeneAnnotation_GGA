#!/bin/bash
usage()
{
    echo "unknown option: ${1}
    usage: fastq_prep_script.sh 
                                 Required inputs:
                                      -r1|--read1 read1.fastq (full path)
                                      -r2|--read2 read2.fastq (full path)
                                      -m|--mintrim mintrim_length (usually 1)
                                      -M|--maxtrim maxtrim_length (usually read length)
                                      -o|--output output_prefix
                                 Options (with defaults shown)
                                      -a1|adapt1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
                                      -a2|adapt2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
                                      
           Script takes paired-end input fastq files (these can be comma-separated lists of fastqs if you want to concatenate multiple fastqs together) and run the lab filter and trim script and then use cutadapt to trim adapter sequences
           
           Right now, this will just work with paired end data, in the future I may add in single end option
           "
}


scriptsPath=
samplename=
mintrim=
maxtrim=

single=0
fileR1=
fileR2=
singleread=

adapt1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapt2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

stranded=
kalStrand=
kalIndex=

while [ "$1" != "" ]; do
    case $1 in
        -r1 | --read1 )      shift
                                fileR1=$1
                                ;;
        -r2 | --read2 )      shift
                                fileR2=$1
                                ;;
        -a1 | --adapt1 )      shift
                                adapt1=$1
                                ;;
        -a2 | --adapt2 )      shift
                                adapt2=$1
                                ;;
        -f | --fastq )      shift
                                singleread=$1
                                ;;
        -o | --outprefix )      shift
                                samplename=$1
                                ;;
        -m | --mintrim )        shift
                                mintrim=$1
                                ;;
        -M | --maxtrim )        shift
                                maxtrim=$1
                                ;;
        -s | --scriptspath )    shift
                                scriptsPath=$1
                                ;;
        --stranded )           shift
                                stranded=$1
                                ;;
        --kalIndex )            shift
                                kalIndex=$1
                                ;;

        --single )    single=1
                                ;;

        -h | --help )           usage "help_message"
                                exit
                                ;;
        * )                     usage $1
                                exit 1
    esac
    shift
done

if [ "$stranded" == "rf-stranded" ]; then
   kalStrand='--rf-stranded'
elif [ "$stranded" == "fr-stranded" ]; then
   kalStrand='--fr-stranded'
elif [ "$stranded" == "unstranded" ]; then
   kalStrand=''
fi



trimandfilterR1=${samplename}.1.trim_${mintrim}_${maxtrim}.minQS_20.fastq
trimandfilterR2=${samplename}.2.trim_${mintrim}_${maxtrim}.minQS_20.fastq
cutadaptR1=${samplename}.1.trim.paired.fastq
cutadaptR2=${samplename}.2.trim.paired.fastq

if [ -n "$singleread" ];then
       echo "single end experiment indicated"
else
       echo "assuming paired-end format"
fi

echo "read1 fastq(s): ${fileR1}"
echo "read2 fastq(s): ${fileR2}"
echo "output prefix: ${samplename}"
echo "trim reads to position ${mintrim} to ${maxtrim}"
echo "Read1 adapter sequence: $adapt1"
echo "Read2 adapter sequence: $adapt2"

echo "Final output filenames: ${samplename}.1.trim.paired.fastq and ${samplename}.2.trim.paired.fastq"


echo "loading modules; gcc/6.2.0; python/2.7.12; cutadapt/1.14"



#### trim and filter
#Check if one or multiple fastqs submitted (per R1 or R2)
   numFastqR1=`echo ${fileR1} | awk 'BEGIN {FS=","} {print NF}'`
   numFastqR2=`echo ${fileR2} | awk 'BEGIN {FS=","} {print NF}'`
   
   #if # of fastqs per read group are uneven, give an error message and terminate the pipeline
   if [ $numFastqR1 -ne $numFastqR2 ]; then
       echo "unequal numbers of fastqs, this is likely a serious problem so terminating the pipeline"
       exit
   fi
   
   #if 1 fastq per read group, then proceed directly to trimand filter. If more than 1, first cat these together, if 0 fastqs then also terminate the pipeline
   if [ "$numFastqR1" -eq 1 ]; then
       echo "one fastq per read group detected - will proceed directly to trimandfilter"
       perl ${scriptsPath}/trim_and_filter_PE.pl -1 ${fileR1} -2 ${fileR2} -a ${mintrim} -b ${maxtrim} -c ${mintrim} -d ${maxtrim} -m 20 -q sanger -o ${samplename}
   elif [ "$numFastqR1" -gt 1 ]; then
       echo "${numFastqR1} fastq per read group detected - will concatenate these and then proceed to trimandfilter" 
       ${scriptsPath}/concatenate_multiple_fastqs_per_sample.sh ${fileR1} ${fileR2} ${samplename}
       perl ${scriptsPath}/trim_and_filter_PE.pl -1 ${samplename}_R1.fastq -2 ${samplename}_R2.fastq -a ${mintrim} -b ${maxtrim} -c ${mintrim} -d ${maxtrim} -m 20 -q sanger -o ${samplename}
   elif [ "$numFastqR1" -eq 0 ]; then
       echo "No fastqs detected - terminating the pipeline"
       exit
   fi



#### Cutadapt trim fastqs
   echo "starting cutadapt trimming"
   cutadapt -f fastq --match-read-wildcards -m 20 -q 10 -a ${adapt1} -A ${adapt2} -o ${cutadaptR1} -p ${cutadaptR2} ${trimandfilterR1} ${trimandfilterR2} &> ../../logs/${samplename}_cutadaptLog.out
   echo "finished cutadapt trimming"
  


## remove the intermediate fastqs
  # rm ${trimandfilterR1}
  # rm ${trimandfilterR2}
#if [ "$numFastqR1" -gt 1 ]; then
  # rm ${samplename}_R1.fastq
  # rm ${samplename}_R2.fastq
#fi

## try adding a sleep step to prevent kallisto problems after this
echo "Kalisto index: ${kalIndex}"
echo "Kalisto strand indo: ${kalStrand}"

echo "Kallisto mapping command to use:  kallisto quant ${kalStrand} -i ${kalIndex} -o $samplename -t 5 -b 30 fastq/${cutadaptR1} fastq/${cutadaptR2}"


echo "starting kallisto mapping"

cd ..


kallisto quant ${kalStrand} -i ${kalIndex} -o $samplename -t 5 -b 30 fastq/${cutadaptR1} fastq/${cutadaptR2}









