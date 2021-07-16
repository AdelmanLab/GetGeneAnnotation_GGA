#!/bin/bash     
            
#Use this if multiple bzipped fastq files exist per sample
#this script takes 3 input parameters
#1: Read1 fastq file paths (comma separated)
#2: Read2 fastq file paths (comma separated)
#3: Sample label for concatenated fastq files 


#split comma-separated fastqs into space-delimited list:
READ1=`echo $1 | sed 's/,/ /g'`
READ2=`echo $2 | sed 's/,/ /g'`
LABEL=$3

#to make this work cleanly first make sure the output files don't already exist (eg. if pipeline fails and then you rerun it again you might cat the fastqs ontop of partial copies from previous piepline run). Here we run some if statements that check if the output file exists and if it does we delete it.
if [ -f ${LABEL}_R1.fastq ] ; then
    rm ${LABEL}_R1.fastq
fi

if [ -f ${LABEL}_R2.fastq ] ; then
    rm ${LABEL}_R2.fastq
fi

#now make the fastq output files
touch ${LABEL}_R1.fastq
touch ${LABEL}_R2.fastq

#loop through R1 fastqs - cat them together in output file. If statements check file extension and use the appropriate cat command.
for FQ_READ1 in $READ1
do
    if [[ "$FQ_READ1" == *".fastq" ]]; then
            cat $FQ_READ1 >> ${LABEL}_R1.fastq
    fi
    
    if [[ "$FQ_READ1" == *".bz2" ]]; then
            bzcat $FQ_READ1 >> ${LABEL}_R1.fastq
    fi
    
    if [[ "$FQ_READ1" == *".gz" ]]; then
            zcat $FQ_READ1 >> ${LABEL}_R1.fastq
    fi
done


#loop through R1 fastqs - cat them together in output file. If statements check file extension and use the appropriate cat command.
for FQ_READ2 in $READ2
do
    if [[ "$FQ_READ2" == *".fastq" ]]; then
            cat $FQ_READ2 >> ${LABEL}_R2.fastq
    fi
    
    if [[ "$FQ_READ2" == *".bz2" ]]; then
            bzcat $FQ_READ2 >> ${LABEL}_R2.fastq
    fi
    
    if [[ "$FQ_READ2" == *".gz" ]]; then
            zcat $FQ_READ2 >> ${LABEL}_R2.fastq
    fi
done



