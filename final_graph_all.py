'''
2023.04.03
Yunseol Park
'''

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import code


def read_file(filename, deam_file, incorrect):
    '''
    File to read in the alignment rates and order them in a graph.
    Args:
        filename (str): name of file that contains different alignment rates (output of find_correctly_aligned.py)
        deam_file (str): name of file that contains the deamination rates
        incorrect (bool/null): a boolean that indicates wheter the alignment is incorrect or correct
                            (default: None - neither incorrect or correct, gives full alignment)
    Returns:
        plot_dict (dict): dictionary that contains different alignment rates and other details needed for plotting
    '''
    # Set the dictionary
    plot_dict = {'Aligner':[], 'Quality':[], 'parameter':[], 'Percentage aligned':[], 'Deamination':[], 'marker':[], 'line':[]}
    deamination = pd.read_csv(deam_file)
    lines = open(filename, 'r').readlines()     # Read alignment rates
    # The alignment rates file contains one line of the file name and one line of different alignment rates
    for l in range(0,len(lines),2):
        aligner, name = lines[l].split('/')[1:]
        aligner = aligner.split('_')[0]     # Get the name of the aligner
        name, q = name.split('.')[:2]       # Get file name and quality score
        if 'deamination' not in name:
            continue
        if 'bwa' in aligner:
            # parameter key refers to mem, aln2, aln1
            if 'mem' in name:
                plot_dict['parameter'].append('-{}'.format(name.split('_')[-1].upper()))
            else:
                plot_dict['parameter'].append(' {}'.format(name.split('_')[-1]))
            aligner = aligner.replace('bwa', 'BWA')
            name = '_'.join(name.split('_')[:-1])
            line = '--'     # Line kwargs for bwa
        else:
            line = '-'
            plot_dict['parameter'].append('')
        if 'alternate' in name:
            aligner += ' ALT'
            name = name.replace('_alternate_allele','')
            mark = 'x'     # marker shape for alternate reads
        else:
            aligner += ' REF'
            mark = '.'
        # After split, use 1 for aligned, 2 for correctly aligned, and 0 for incorrectly aligned
        if incorrect:
            data = float(lines[l+1].split()[0])
        elif incorrect == False:
            data = float(lines[l+1].split()[2])
        else:
            data = float(lines[l+1].split()[1])
        # Deamination level
        level = float(deamination[deamination['Sample']==name.split('_')[2]]['C to T'])
        # Fill out the rest of the information
        plot_dict['Aligner'].append(aligner)
        plot_dict['Quality'].append(q)
        plot_dict['Percentage aligned'].append(data)
        plot_dict['Deamination'].append(level)
        plot_dict['marker'].append(mark)
        plot_dict['line'].append(line)
    return plot_dict

def plot_full(filename, deam_file='data/deamination.csv', incorrect=False):
    '''
    Args:
        filename (str): name of file that contains different alignment rates (output of find_correctly_aligned.py)
        deam_file (str): name of file that contains the deamination rates
        incorrect (bool/null): a boolean that indicates wheter the alignment is incorrect or correct
                            (default: None - neither incorrect or correct, gives full alignment)
    Returns:
        None
    '''
    plot_dict = read_file(filename, deam_file, incorrect)
    plot_df = pd.DataFrame(plot_dict)
    # Create subplots and set parameters
    figure, axes = plt.subplots(4, 3, figsize=(30,30))
    figure.tight_layout(pad=12)
    plt.rcParams.update({'font.size': 35, "font.family": "serif", "mathtext.default": "regular", "mathtext.fontset": "dejavuserif"})
    # Y-axis limit
    if incorrect:
        lim = [0, 0.5]
        subtitle = [7.5, 0.64]
        axis = 'incorrectly '
    else:
        lim = [80, 100]
        subtitle = [7, 101.45]
        if incorrect == False:
            axis = 'correctly '
        else:
            axis = ''
    # Axis labels
    figure.text(0.5, 0.02, 'Deamination (%)', ha='center', va='center', fontsize=40)
    figure.text(0.015, 0.5, 'Percentage {}aligned (%)'.format(axis), ha='center', va='center', rotation='vertical', fontsize=40)
    # counters for axes
    row = 0
    col = 0
    # Colors of graph in HEX (order: vg REF, vg ALT, bwa REF, bwa ALT)
    #'#b2182b', '#ef8a62', '#2166ac', '#67a9cf' ['#762a83', '#af8dc3', '#1b7837', '#7fbf7b']
    colors = ['#af8dc3', '#762a83', '#7fbf7b', '#1b7837']
    # Label for the subplots and counter
    label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)']
    counter = 0

    for vg_quality in ['q50', 'q60']:
        for bwa_quality in ['q25', 'q30']:
            for parameter in [' aln1', ' aln2', '-MEM']:
                if parameter == '-MEM':
                    if bwa_quality == 'q25':
                        bwa_quality = 'q50'
                    else:
                        bwa_quality = 'q60'
                needed = ['vg REF', 'vg ALT', 'BWA REF', 'BWA ALT']
                quality = [vg_quality, vg_quality, bwa_quality, bwa_quality]
                p_full = ['', '', parameter, parameter]
                for i, (sub, q, p) in enumerate(zip(needed, quality, p_full)):
                    # Subset only one type of alignment
                    sub_df = plot_df.loc[(plot_df['Aligner']==sub) & (plot_df['Quality']==q) & (plot_df['parameter']==p)]
                    # Plot graph
                    sns.regplot(ax=axes[row,col], x='Deamination', y='Percentage aligned', data=sub_df, color=colors[i],
                               marker=sub_df['marker'].iloc[0], line_kws={'ls':sub_df['line'].iloc[0]}, ci=None, label=sub, scatter_kws={'s':150})
                # Plot subtitles and limits
                axes[row,col].text(subtitle[0], subtitle[1],"{}    vg  {} vs BWA{} {}".format(label[counter], vg_quality, parameter, bwa_quality))
                axes[row,col].set_ylim(lim)
                if incorrect:
                    axes[row,col].legend(fontsize=25, loc='upper left')
                else:
                    axes[row,col].legend(fontsize=25, loc='lower left')
                axes[row,col].xaxis.set_tick_params(labelsize=25)
                axes[row,col].yaxis.set_tick_params(labelsize=25)
                axes[row,col].set_xlabel(None)
                axes[row,col].set_ylabel(None)
                col += 1
                if col > 2:
                    col = 0
                    row += 1
                counter += 1
    plt.savefig('images/comparison_correct.png')

plot_full('aligned_percentage_final_wtolerance2.txt')