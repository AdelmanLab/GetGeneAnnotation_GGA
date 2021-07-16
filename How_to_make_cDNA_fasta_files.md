
## Making cDNA fasta files from gtfs

Here is example code, starting from an ENSEMBL GTF, to filter for basic annotations and then to make the gene-info text files and cDNA fasta files. Adjust your starting GTF and filenames as appropiate.

Note that the code to make the transcript to genes info text file should work for ENSEMBL annotations, but the order of attributes in the 9th column differs for other GTF sources (eg. Genecode) and so the code will have to be adjusted for other GTFs from other sources.

```
# starting with a gzipped gtf file: Mus_musculus.GRCm38.99.ntc.filtered.gtf.gz

#Filter for basic transcripts only:
zcat Mus_musculus.GRCm38.99.ntc.filtered.gtf.gz | grep 'tag "basic"' > temp.mm10.basic.gtf
zcat Mus_musculus.GRCm38.99.ntc.filtered.gtf.gz | awk '{if($3=="gene"){print $0}}' | cat - temp.mm10.basic.gtf | sort -k 1,1 -k 4,4n > mm10.basic.gtf
rm temp.mm10.basic.gtf

#check how many transcripts in original vs filtered GTF
zcat Mus_musculus.GRCm38.99.ntc.filtered.gtf.gz | awk '{if($3=="transcript"){print $0}}' | wc -l
awk '{if($3=="transcript"){print $0}}' mm10.basic.gtf | wc -l
#This reduces transcripts from 142,100 to 81,371

#Pull out a text file with transcriptID, coordinates, and gene info (ID, name, and biotype)
#apologies for the gem of a one-liner here...
cat mm10.basic.gtf | sed '1,5d' | awk '{if($3=="transcript"){print $0}}' | sed 's/gene_id /;/g' | sed 's/transcript_id /;/g' | sed 's/gene_name /;/g' | sed 's/gene_biotype /;/g' | awk 'BEGIN {FS=";"};{OFS="\t";print($1,$5,$2,$8,$11)}' | sed 's/"//g' | awk '{OFS="\t";print($9,$1,$4,$5,$7,$10,$11,$12,$13)}' | sort -k 1,1 > mm10.basic.transcripts.to.genes.info.txt

#Get cDNA for basic transcripts
#download gffread tools
    wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.11.8.Linux_x86_64.tar.gz
#Unzip gffread tools
    tar -xzf gffread-0.11.8.Linux_x86_64.tar.gz
#extract DNA sequences for spliced exons using the lab mm10 fasta file (dowloaded from the UCSC browser)
    gffread-0.11.8.Linux_x86_64/gffread -w mm10.basic.cDNA.fa -g mm10.fa mm10.basic.gtf
#Remove gffread tools files    
    rm -r gffread-0.11.8.Linux_x86_64
    rm gffread-0.11.8.Linux_x86_64.tar.gz

#Make cDNA kallisto index
   module load gcc/6.2.0
   module load kallisto/0.45.1
   kallisto index -i mm10.basic.cDNA.idx mm10.basic.cDNA.fa   

```


