#!/bin/bash
#SBATCH -p batch
###SBATCH -p express
###SBATCH -t 60:00:00
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=gabriela.toomer@okstate.edu
#SBATCH --mail-type=end
module load bowtie2/2.3.4.1 fastqc/0.11.7 samtools/1.10 trimmomatic/0.38 cufflinks/2.2.1 tophat/2.1.1

#!/bin/bash
# Basic for loop
names= 'T0_1 T0_2 T0_3
        M'

for name in $names
do

tophat2 --num-threads 32 -o /scratch/gatoo/output_"$name"  -G /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf
  --transcriptome-index /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/
  --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1
  --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000
  /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Sequence/Bowtie2Index/genome
  /scratch/gatoo/"$name"_F.fastq
  /scratch/gatoo/"$name"_R.fastq

cufflinks --num-threads 32 -o /scratch/gatoo/output_cuff_"$name"
  --frag-bias-correct /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Sequence/WholeGenomeFasta/genome.fa
  --multi-read-correct -G /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf
  /scratch/gatoo/output_"$name"/accepted_hits.bam

  echo /scratch/gatoo/output_"$name"_cuff/genes.fpkm_tracking
  done
  echo All done


#gunzip /scratch/gatoo/T2_2_*.fastq.gz
#fastqc *.fastq.gz

##Build the genome annotation
#bowtie2-build dna_index/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz Bos_taurus
#tophat2 --num-threads 32 -o /scratch/gatoo/T2_2 BHV1 /scratch/gatoo/T2_2_F.fastq /scratch/gatoo/T2_2_R.fastq
#cufflinks BHV1 /scratch/gatoo/T2_1_b
#Get counts using cufflinks
#cufflinks --num-threads 32 -o T2_2 --frag-bias-correct BHV1.fa --multi-read-correct -G BHV1.gtf /scratch/gatoo/T2_2/accepted_hits.bam

##After downloading from igenomes, Ready-To-Use Reference Sequences and Annotations
#http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/
#For Caw http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Bos_taurus/UCSC/bosTau8/Bos_taurus_UCSC_bosTau8.tar.gz
#For Mouse http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
