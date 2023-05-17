#!/bin/bash

# Write a log file to record the time, memory, cpu usage, and errors
exec 1> bwa_log_final.out 2>&1

# Path to store the alignment results
RESULT=data/bwa_aligned_final
# Path and name of the reference genome
REF=data/reference/hs37d5.11.fa
# Path to the input data (fq.gz)
PREFIX=data/simulated_data

# If the result directory does not exist, make directory
if [ ! -d $RESULT ]
then
    mkdir $RESULT
fi
# If the reference has not been indexed, run bwa index
if [ ! -e $REF.amb ]
then
    printf "creating index for ${REF##data*/}\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bwa/bwa index $REF
fi

for fq_path in $PREFIX/*
do
    # Set names for each alignment
    echo $fq_path
    NAME=${fq_path##data*50bp_}
    NAME=${NAME%.fq.gz}
    # Code below can be used to remove certain files from alignment
    #if [ $NAME == sliding_window ] || [ $NAME == sliding_window_alternate_allele ] || [[ $NAME == *allele_RISE00* ]]
    #then
    #    continue
    #fi

    # Run alignment for bwa aln with -n 0.02 parameter
    printf "start bwa aln 0.02 of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bwa/bwa aln -l 1024 -n 0.02 -q 15 $REF $fq_path > $RESULT/$NAME\_aln2.sai
    printf "finished aln 0.02\n"
    # bwa samse converts sai to sam file format
    printf "start bwa aln 0.02 samse of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bwa/bwa samse $REF $RESULT/$NAME\_aln2.sai $fq_path> $RESULT/$NAME\_aln2.sam
    printf "finished aln 0.02\n"
    # Convert sam to bam file format
    printf "start bwa aln 0.02 samtools view of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bin/samtools view -b $RESULT/$NAME\_aln2.sam > $RESULT/$NAME\_aln2.bam
    printf "finished aln 0.02\n"
    # Sort the file for samtools markdup
    bin/samtools sort -o $RESULT/$NAME\_aln2_sorted.bam $RESULT/$NAME\_aln2.bam
    # Remove PCR duplicates
    bin/samtools markdup -r -s $RESULT/$NAME\_aln2_sorted.bam $RESULT/$NAME\_aln2_rmdup.bam
    #cd $current_dir
    # Index files are required for samtools view with -q parameter (quality score)
    printf "make index"
    bin/samtools index $RESULT/$NAME\_aln2_rmdup.bam
    printf "index made"
    # select alignments above quality score 25 and 30
    bin/samtools view -b -q 25 $RESULT/$NAME\_aln2_rmdup.bam > $RESULT/$NAME\_aln2_q25.bam
    bin/samtools view -b -q 30 $RESULT/$NAME\_aln2_rmdup.bam > $RESULT/$NAME\_aln2_q30.bam
    # Index the final file for visualization
    bin/samtools index -M $RESULT/$NAME\_aln2_q30.bam
    
    # Run alignment for bwa aln with -n 0.01 -o 2 parameter
    printf "start bwa aln 0.01 of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bwa/bwa aln -l 1024 -n 0.01 -o 2 -q 15 $REF $fq_path > $RESULT/$NAME\_aln1.sai
    printf "finished aln 0.01\n"
    # Convert sai to sam
    printf "start bwa aln 0.01 samse of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bwa/bwa samse $REF $RESULT/$NAME\_aln1.sai $fq_path > $RESULT/$NAME\_aln1.sam
    printf "finished aln 0.01\n"
    # Convert sam to bam
    printf "start bwa aln 0.01 samtools view of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bin/samtools view -b $RESULT/$NAME\_aln1.sam > $RESULT/$NAME\_aln1.bam
    printf "finished aln 0.01\n"
    # Sort and remove PCR duplicates
    bin/samtools sort -o $RESULT/$NAME\_aln1_sorted.bam $RESULT/$NAME\_aln1.bam
    bin/samtools markdup -r -s $RESULT/$NAME\_aln1_sorted.bam $RESULT/$NAME\_aln1_rmdup.bam
    # Make index and select for quality scores above 25 and 30
    printf "make index"
    bin/samtools index $RESULT/$NAME\_aln1_rmdup.bam
    printf "index made"
    bin/samtools view -b -q 25 $RESULT/$NAME\_aln1_rmdup.bam > $RESULT/$NAME\_aln1_q25.bam
    bin/samtools view -b -q 30 $RESULT/$NAME\_aln1_rmdup.bam > $RESULT/$NAME\_aln1_q30.bam
    # Index the final file for visualization
    bin/samtools index -M $RESULT/$NAME\_aln1_q30.bam

    # Run alignment for bwa mem with default parameters
    printf "start bwa mem of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bwa/bwa mem $REF $fq_path > $RESULT/$NAME\_mem.sam
    printf "finished mem\n"
    # Convert sam to bam
    printf "start bwa mem samtools view of $NAME\n"
    /usr/bin/time -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' bin/samtools view -b $RESULT/$NAME\_mem.sam > $RESULT/$NAME\_mem.bam
    printf "finished mem\n"
    # Sort and remove PCR duplicates
    bin/samtools sort -o $RESULT/$NAME\_mem_sorted.bam $RESULT/$NAME\_mem.bam
    bin/samtools markdup -r -s $RESULT/$NAME\_mem_sorted.bam $RESULT/$NAME\_mem_rmdup.bam
    # Make index and select for those with quality scores above 25 and 30
    printf "make index"
    bin/samtools index $RESULT/$NAME\_mem_rmdup.bam
    printf "index made"
    bin/samtools view -b -q 50 $RESULT/$NAME\_mem_rmdup.bam > $RESULT/$NAME\_mem_q25.bam
    bin/samtools view -b -q 60 $RESULT/$NAME\_mem_rmdup.bam > $RESULT/$NAME\_mem_q30.bam
    # Index the final file for visualization
    bin/samtools index -M $RESULT/$NAME\_mem_q30.bam
done

