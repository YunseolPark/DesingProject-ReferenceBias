#!/bin/bash

# Write a log file to record the time, memory, cpu usage, and errors
exec 1> vg_script_final.out 2>&1

# Path to store the alignment results
RESULT=data/vg_aligned_final
# Path and name of the reference genome
REF=data/reference/hs37d5.11.fa
VAR=data/reference/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
# Index path
INDEX=data/reference/final
# Path to the input data (fq)
PREFIX=data/simulated_data

# If the result directory does not exist, make directory
if [ ! -d $RESULT ]
then
    mkdir $RESULT
fi

# If the graph has not been constructed, construct graph
printf "creating graph\n"
#/usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' /opt/vg construct -C -R 11 -r $REF -v $VAR -m 32 > $INDEX.vg
/usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' /opt/vg construct -r $REF -v $VAR -m 32 > $INDEX.vg

# If the reference has not been indexed, make xg index
printf "creating xg index\n"
/usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' /opt/vg index -x $INDEX.xg $INDEX.vg

# If the reference has not been indexed, make gcsa index
printf "creating gcsa index\n"
/usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' /opt/vg prune -r $INDEX.vg > $INDEX.pruned.vg
/usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' /opt/vg index -g $INDEX.gcsa $INDEX.pruned.vg

for fq_path in $PREFIX/*
do
    # Set names for each alignment
    echo $fq_path
    NAME=${fq_path##data*50bp_}
    NAME=${NAME%.fq.gz}

    printf "Map $NAME"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' /opt/vg map -f $fq_path -x $INDEX.xg -g $INDEX.gcsa --surject-to bam -k 15 -w 1024 > $RESULT/$NAME.bam
    # Sort the file for samtools markdup
    bin/samtools sort -o $RESULT/$NAME.sorted.bam $RESULT/$NAME.bam
    # Remove PCR duplicates
    bin/samtools markdup -r -s $RESULT/$NAME.sorted.bam $RESULT/$NAME.rmdup.bam
    #cd $current_dir
    # Index files are required for samtools view with -q parameter (quality score)
    bin/samtools index $RESULT/$NAME.rmdup.bam
    # select alignments above quality score 25 and 30
    bin/samtools view -b -q 50 $RESULT/$NAME.rmdup.bam > $RESULT/$NAME.q50.bam
    bin/samtools view -b -q 60 $RESULT/$NAME.rmdup.bam > $RESULT/$NAME.q60.bam
    # Save sam file from the resultant bam file to check for differences
    #bin/samtools view -h $RESULT/$NAME.q50.bam > $RESULT/$NAME.q50.sam
    #bin/samtools view -h $RESULT/$NAME.q60.bam > $RESULT/$NAME.q60.sam
    # Index the final file for visualization
    #bin/samtools index -M $RESULT/$NAME.q60.bam
done