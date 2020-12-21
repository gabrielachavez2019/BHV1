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
#module load tophat/2.1.1

Lat_2_F.fastq.gz

#gunzip /scratch/gatoo/T2_2_*.fastq.gz
#fastqc *.fastq.gz

bowtie2-build genome/Bos_taurus_UCSC_bosTau8.tar.gz /scracth/gatoo/genome/Bos_taurus
#bowtie2-build /scratch/gatoo/dna/Bos_taurus.ARS-UCD1.2.dna.primary_assembly.*.fa.gz /scracth/gatoo/Bos_taurus
#Bos_taurus.ARS-UCD1.2.101.gff3.gz
#tophat2 --num-threads 32 -o /scratch/gatoo/T2_2 BHV1 /scratch/gatoo/T2_2_F.fastq /scratch/gatoo/T2_2_R.fastq
#cufflinks BHV1 /scratch/gatoo/T2_1_b
#Get counts using cufflinks
#cufflinks --num-threads 32 -o T2_2 --frag-bias-correct BHV1.fa --multi-read-correct -G BHV1.gtf /scratch/gatoo/T2_2/accepted_hits.bam

#tophat2 -o output -G Homo_sapiens/GFF/ref_GRCh38.p12_top_level.gff3 --transcriptome-index Homo_sapiens/RNA/ --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1 --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000 ../rna_hs19 ../trimmed_T0-a.fastq
#tophat2 -o output -G /scratch/gatoo/bos_taurus/Bos_taurus.ARS-UCD1.2.101.gff3.gz --transcriptome-index Bos_taurus/cdna/ --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1 --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000 /scratch/gatoo/bos_taurus/cdna  /scratch/gatoo/[i].fastq
#bowtie2-build --num-threads 32 /scratch/gatoo/bos_taurus/dna_index/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz /scratch/gatoo/Bos_taurus/Bos_taurus

bowtie2-build dna_index/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz Bos_taurus


####
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
tophat2  --num-threads 32  -o /scratch/gatoo/output_v_m_1 -G BHV1.gtf BHV1 /scratch/gatoo/m_1_F.fastq.gz  /scratch/gatoo/m_1_R.fastq.gz

#!/bin/bash
# Basic for loop
names= 'T0_1 T0_2 T0_3
        T1_1 T1_2
        T2_1 T2_2 T2_3
        T3_1 T3_2 T3_3
        M'

for name in $names
do
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
tophat2 --num-threads 32 -o /scratch/gatoo/output  -G /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf --transcriptome-index /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/  --read-edit-dist 2 --read-gap-length 2 --read-mismatches 2 --read-realign-edit-dist 1 --library-type fr-unstranded --min-intron-length 50 --max-intron-length 300000 /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Sequence/Bowtie2Index/genome  /scratch/gatoo/T2_1_F.fastq /scratch/gatoo/T2_1_R.fastq
cufflinks --num-threads 32 -o /scratch/gatoo/output --frag-bias-correct /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Sequence/WholeGenomeFasta/genome.fa --multi-read-correct -G /scratch/gatoo/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf /scratch/gatoo/output/accepted_hits.bam


cufflinks -o /scratch/gatoo/HSV1-mouse/output_cuff_HSV1_Fem-WT/ --frag-bias-correct /scratch/gatoo/HSV1-mouse/HSV1.fasta --multi-read-correct -G /scratch/gatoo/HSV1-mouse/HSV1.gff3 /scratch/gatoo/HSV1-mouse/output_HSV1_Fem-WT/accepted_hits.bam
