#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=gabriela.toomer@okstate.edu
#SBATCH --mail-type=end
module load bowtie2/2.3.4.1 fastqc/0.11.7 samtools/1.10 trimmomatic/0.38 cufflinks/2.2.1 tophat/2.1.1 R

##Concatenate files Mocks or uninfected:
  echo Begining concatenating files...

cat /scratch/gatoo/SecNGS/LIB109080_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-dLAT_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109080_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-dLAT_R.fastq.gz

cat /scratch/gatoo/SecNGS/LIB109081_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-WT_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109081_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-WT_R.fastq.gz

cat /scratch/gatoo/SecNGS/LIB109082_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-U1_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109082_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Fem-U1_R.fastq.gz

cat /scratch/gatoo/SecNGS/LIB109083_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109083_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_R.fastq.gz

cat /scratch/gatoo/SecNGS/LIB109084_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-dLAT_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109084_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-dLAT_R.fastq.gz

cat /scratch/gatoo/SecNGS/LIB109085_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-WT_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109085_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-WT_R.fastq.gz

#BoHSV1-cow

echo Finishing full size library files...

names='Fem-dLAT Fem-WT Fem-U1
       Male-dLAT Male-WT Male-U1'

for name in $names
 do
  echo Starting mapping of "$name"
tophat2 --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_"$name"  -G /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf --transcriptome-index /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/ --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1 --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000 /scratch/gatoo/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome /scratch/gatoo/HSV1-mouse/"$name"_F.fastq.gz /scratch/gatoo/HSV1-mouse/"$name"_R.fastq.gz

cufflinks --num-threads 32 -o /scratch/gatoo/HSV1-mouse/output_cuff_"$name"  --frag-bias-correct /scratch/gatoo/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa  --multi-read-correct -G /scratch/gatoo/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf /scratch/gatoo/HSV1-mouse/output_"$name"/accepted_hits.bam

cp  /scratch/gatoo/HSV1-mouse/output_cuff_"$name"/genes.fpkm_tracking HSV1-mouse/"$name".txt

###Assemlbies file should have separate lines with the full path to the gtf files per sample
/scratch/gatoo/output_cuff_M/transcripts.gtf
/scratch/gatoo/output_cuff_T0_2/transcripts.gtf
/scratch/gatoo/output_cuff_T2_2/transcripts.gtf
/scratch/gatoo/output_cuff_T0_1/transcripts.gtf
/scratch/gatoo/output_cuff_T0_3/transcripts.gtf
##Then run
  done
echo All done

cat /scratch/gatoo/SecNGS/LIB109083_S*_R1_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_F.fastq.gz
cat /scratch/gatoo/SecNGS/LIB109083_S*_R2_001.fastq.gz >> /scratch/gatoo/HSV1-mouse/Male-U1_R.fastq.gz


names='m t0'
       #t1 t2 t3'
for name in $names
 do
  echo

cuffmerge -g /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf -s /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.fa -p 8 assemblies.txt
#Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate
cuffdiff -o diff_out -b /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.fa -p 8 -L m,t0 -u merged_asm/merged.gtf \/scratch/gatoo/output_m1/accepted_hits.bam,/scratch/gatoo/output_m_2/accepted_hits.bam,/scratch/gatoo/output_m_3/accepted_hits.bam   \/scratch/gatoo/output_t0_2/accepted_hits.bam,/scratch/gatoo/output_t0_1/accepted_hits.bam,/scratch/gatoo/output_t0_3/accepted_hits.bam
#Create the plots with the R script
Rscript Plots.R

  echo "$name" counting: done!
 done
echo All done
