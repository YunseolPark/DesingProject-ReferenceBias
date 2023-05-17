'''
2023.04.03
Yunseol Park
'''

import sys
import os
import subprocess
import code

def calculate_aln(aln_file, write_file, tolerance):
    '''
    Function to calculate the different alignment rates.
    Args:
        aln_file (str): name of the BAM file
        write_file (TextIOWrapper): the opened file to write details in
        tolerance (int): the range of tolerance (default: None)
    Returns:
        None
    '''
    write_name = aln_file
    if '_q' in aln_file:
        write_name = aln_file.replace('_q', '.q')
    if 'mem' in aln_file:
        if '30' in aln_file:
            write_name = aln_file.replace('q30', 'q60')
        else:
            write_name = aln_file.replace('q25', 'q50')
    print(write_name, file=write_file)
    # Read in the aligned files usign samtools -> list of each line
    content = subprocess.check_output(['bin/samtools','view', aln_file])
    content = content.decode("utf-8").split('\n')
    correct = {}
    unique_aln = set()
    incorrect = {}

    for line in content:
        # Last line is an empty string
        if line == '':
            break
        # Find the position of the read (reference) - 1st column in the sam file
        ref_pos = int(line.split('\t')[0][6:].split('-')[0])
        # Find the position of the mapping/alignment - 4th column in the sam file
        mapped_pos = int(line.split('\t')[3])
        # Add the positions to the dictionary if the two positions are the same
        unique_aln.add(ref_pos)
        if tolerance:
            check = ref_pos - tolerance <= mapped_pos <= ref_pos + tolerance
        else:
            check = ref_pos == mapped_pos
        # If there is already a correct alignment at that position, keep the correct alignment.
        if ref_pos in correct:
            correct[ref_pos].append(line)
        elif check:
            correct[ref_pos] = [line]
        else:   # Incorrect alignment
            if ref_pos in incorrect:
                incorrect[ref_pos].append(line)
            else:
                incorrect[ref_pos] = [line]

    correctly_aln = len(correct.keys())     # Correctly aligned reads
    aln = len(unique_aln)       # Aligned reads

    # 1502550: total number of reads in the simulated data
    total = 1502550
    percentage_aln = aln/total * 100      # Percentage of aligned reads (out of the entire reads)
    percentage_non = (len(incorrect.keys()))/total * 100    # Percentage of reads that are not aligned (out of the entire reads)
    percentage_correct = correctly_aln/total * 100        # Percentage of the correctly aligned reads (out of the entire reads)
    correct_outof_aln = correctly_aln/aln * 100     # Percentage of correctly aligned reads from the aligned reads

    # Write to file
    print(percentage_non, percentage_aln, percentage_correct, correct_outof_aln, file=write_file)


def filter_exact_match(sim_dir, output, vg_dir=None, bwa_dir=None, tolerance=None):
    '''
    Args:
        sim_dir (str): directory name of the simulated data
        output (str): name of file to write the alignment rates to
        vg_dir (str): directory name for the vg alignment files
        bwa_dir (str): directory name for the BWA alignment files
        tolerance (int): the range of tolerance (default: None)
    Returns:
        None
    '''
    write_file = open(output, 'w')
    for filename in os.listdir(sim_dir):
        #print(filename)
        if '.fq.gz' not in filename:
            continue
        # Get the names of the aligned files and calculate the alignment percentages
        aligned_file = filename.replace('HO_chr11_50bp_', '')
        if vg_dir:
            vg60_file = os.path.join(vg_dir, aligned_file.replace('fq.gz', 'q60.bam'))
            calculate_aln(vg60_file, write_file, tolerance)
            vg50_file = os.path.join(vg_dir, aligned_file.replace('fq.gz', 'q50.bam'))
            calculate_aln(vg50_file, write_file, tolerance)
        if bwa_dir:
            aln1_30_file = os.path.join(bwa_dir, aligned_file.replace('.fq.gz', '_aln1_q30.bam'))
            calculate_aln(aln1_30_file, write_file, tolerance)
            aln1_25_file = os.path.join(bwa_dir, aligned_file.replace('.fq.gz', '_aln1_q25.bam'))
            calculate_aln(aln1_25_file, write_file, tolerance)
            aln2_30_file = os.path.join(bwa_dir, aligned_file.replace('.fq.gz', '_aln2_q30.bam'))
            calculate_aln(aln2_30_file, write_file, tolerance)
            aln2_25_file = os.path.join(bwa_dir, aligned_file.replace('.fq.gz', '_aln2_q25.bam'))
            calculate_aln(aln2_25_file, write_file, tolerance)
            mem60_file = os.path.join(bwa_dir, aligned_file.replace('.fq.gz', '_mem_q30.bam'))
            calculate_aln(mem60_file, write_file, tolerance)
            mem50_file = os.path.join(bwa_dir, aligned_file.replace('.fq.gz', '_mem_q25.bam'))
            calculate_aln(mem50_file, write_file, tolerance)
    write_file.close()
        
filter_exact_match('data/simulated_data', 'aligned_percentage_final_wtolerance2.txt', 'data/vg_aligned_final', 'data/bwa_aligned_final', 5)
