#!/bin/bash
# Basic for loop
names= 'T0_1 T0_2 T0_3
        T1_1 T1_2
        T2_1 T2_2 T2_3
        T3_1 T3_2 T3_3
        M'

for name in $names
do
echo   /scratch/gatoo/output_cuff_"$name"/transcripts.gtf /scratch/gatoo/output_cuff_"$name"/transcripts.gtf
  /scratch/gatoo/output_cuff_"$name"/transcripts.gtf /scratch/gatoo/output_cuff_"$name"/transcripts.gtf

done
echo All done
