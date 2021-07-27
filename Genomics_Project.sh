
#Check the quality of the data using FASTQC--------------------------------------------------------------

conda activate    #All tools were installed under miniconda 3 in the base environment

fastqc SRR576933.fastq    #control

fastqc SRR576938.fastq    #Experimental condition


#Mapping-------------------------------------------------------------------------------------------------

bowtie-build EC_K12.fna E_coli_K12    #Building the genome index file to perform the mapping later

#Mapping the ChIp readings with the reference genome using the genome index file previously built

bowtie E_coli_K12 -q SRR576933.fastq  -v 2 -m 1 -3 1 -S 2> SRR576933.out > SRR576933.sam 

bowtie E_coli_K12 -q SRR576938.fastq  -v 2 -m 1 -3 1 -S 2> SRR576938.out > SRR576938.sam


#Verify Two ENCODE quality metrics-----------------------------------------------------------------------

#For the experiment

samtools view -bS SRR576933.sam  > SRR576933.bam

samtools view -F 0x0204 -o - SRR576933.bam | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > SRR576933_experiment.tagAlign.gz

Rscript ./run_spp.R -c=SRR576933_experiment.tagAlign.gz  -savp -out=SRR576933_experiment_phantompeaks

#for the control

samtools view -bS SRR576938.sam  > SRR576938.bam

samtools view -F 0x0204 -o - SRR576938.bam | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > SRR576938_experiment.tagAlign.gz

Rscript ./run_spp.R -c=SRR576938_experiment.tagAlign.gz  -savp -out=SRR576938_experiment_phantompeaks


#Peak calling using MACS2-------------------------------------------------------------------------------

macs2 callpeak -t SRR576933.sam -c SRR576938.sam --format SAM  --gsize 4641652 --name "macs2"  --bw 400 --keep-dup 2 --bdg --nomodel --extsize 200 &> MACS.out


#Motif analysis-----------------------------------------------------------------------------------------

# bedtools getfasta -fi EC_K12.fna -bed macs2_summits.bed -fo macs2_peaks.fa  #from the peaks only

bedtools getfasta -fi EC_K12.fna -bed macs2_peaks.narrowPeak -fo macs2_NarrowPeaks.fa #from the bed+4 file
 
perl -lane '$start=$F[1]-100 ; $end = $F[2]+100 ; print "$F[0]\t$start\t$end"' macs2_summits.bed > macs2_summits+-100.bed #from the summits +-100 bp

bedtools getfasta -fi EC_K12.fna -bed macs2_summits+-100.bed -fo macs2_summits+-100.fa


