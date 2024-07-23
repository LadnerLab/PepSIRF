#-------DELETE LATER--------
import sys
sys.path.append("/Users/scg283/Documents/GitHub/Modules")
#---------------------------

import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
import fastatools as ft
from scipy.signal import find_peaks

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i', '--input-dir',  help='Directory with alignment output files and files that contain the mapped location of peptides', required=True)
    parser.add_argument('--window-size', type=int, default=30,  help='Size of AA window to use for identifying core epitopes.', required=False)
    parser.add_argument('--max-zeros', type=int, default=5, help='Maximum number of zero counts a window can contain.', required=False)
    parser.add_argument('--max-overlap', type=int, default=8, help='Maximum AA overlap a window can have with a previously selected window.', required=False)
    parser.add_argument('--peptide-overlap', type=int, default=9, help='Peptide sequence should overlap at least this amount to be included in output data.', required=False)
    parser.add_argument('--peak-overlap-window-size', type=int, default=10, help='Window size around found peaks in which containing peptides will be removed for the next iteration.', required=False)
    parser.add_argument('--include-iter-vis', action="store_true", help='Output each chart given the iteration.', required=False)
    parser.add_argument('-o', '--output-dir', default="find_epitopes_out", help='Name of directory to output line plots.')
    
    args = parser.parse_args()

    directory_path = args.input_dir
    alignment_to_use_dict = read_check_align_file(directory_path)

    # get original align counts and peptide positions 
    alignCountsD, file_2_pep_pos_dict = process_files_probes(
                                            probes_dict=alignment_to_use_dict, 
                                            directory_path=directory_path
                                            )

    #windows = find_core_epitopes(alignCountsD, args.window_size, args.max_zeros, args.max_overlap)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    else:
        print(f"Warning: The directory \"{args.output_dir}\" already exists. Files may be overwritten.\n")

    windows = iterative_peptide_finder(alignment_to_use_dict, directory_path, args.window_size, args.max_zeros, args.max_overlap, 
                                                                args.peak_overlap_window_size, args.output_dir, args.include_iter_vis)

    generate_out_data(args.output_dir, directory_path, windows, file_2_pep_pos_dict, args.peptide_overlap)

    create_line_charts(alignCountsD, windows, args.output_dir)


def iterative_peptide_finder(alignment_to_use_dict, directory_path, window_size, max_zeros, max_overlap, 
                                                                    peak_ovlp_window_size, output_dir, include_iter_vis):
    out_dict = dict()

    for filename, data in alignment_to_use_dict.items():
        if include_iter_vis:
            iter_vis_out = os.path.join(output_dir, f"{filename.split('_')[-2]}_iterations")
            if not os.path.exists(iter_vis_out):
                os.mkdir(iter_vis_out)

        iter_num = 1
        windows = list()
        removed_peptides = set()

        aligned_probes_path = os.path.join(directory_path, filename.replace('.fasta', '_probesAligned.txt'))

        peak_found = True
        window_found = True
        while peak_found:
            peak_found = False
            score = 0

            # generate scores without used peptides
            alignCountsD, pep_pos_dict = process_file_probes(
                                                    data=data, 
                                                    aligned_probes_path=aligned_probes_path,
                                                    removed_peptides=removed_peptides
                                                    )

            y = list(alignCountsD.values())
            
            # find highest peak
            peaks, _ = find_peaks(y)
            peaks = list(peaks)

            if len(peaks) > 0:
                valid_window = False
                peaks_empty = False
                while not valid_window and not peaks_empty:
                    # get max peak (center if there are multiple)
                    max_peak = peaks[np.argmax([y[peak] for peak in peaks])]
                    max_peaks = list()
                    max_peak_indices = list()
                    for peak_idx in range(len(peaks)):
                        peak = peaks[peak_idx]
                        if y[peak] == y[max_peak]:
                            max_peaks.append(peak)
                            max_peak_indices.append(peak_idx)
                    max_peak_idx = len(max_peaks)//2
                    max_peak = max_peaks[max_peak_idx]

                    # get designed peptide window
                    left_border, right_border = generate_window(max_peak, window_size, int(data))

                    # test if window passes thresholds
                    if y[left_border:right_border].count(0) <= max_zeros and \
                                    all(get_overlap((left_border, right_border), (x[0], x[1])) <= max_overlap for x in windows):
                        valid_window = True
                    else:
                        # remove possible mass peak
                        peaks.pop(max_peak_indices[max_peak_idx])
                        if len(peaks) == 0:
                            peaks_empty = True

                if not peaks_empty:
                    peak_found = True

                    # get the overlap window
                    pep_ovlp_win = (max_peak - (peak_ovlp_window_size // 2), max_peak + (peak_ovlp_window_size // 2))

                    # remove overlapping peptides
                    for pep, pep_coords in pep_pos_dict.items():
                        if get_overlap(pep_ovlp_win, pep_coords) > 0:
                            removed_peptides.add(pep)

                    windows.append((left_border, right_border))

                    if include_iter_vis:
                        create_line_chart(
                                x = list(alignCountsD.keys()),
                                y = list(alignCountsD.values()),
                                windows = windows,
                                title = filename,
                                out_dir = os.path.join(iter_vis_out, f"step_{iter_num}.png")
                                )

                    iter_num += 1

        out_dict[filename] = windows

    return out_dict


def generate_window(max_peak, window_size, max_x):
    left_border = max_peak - (window_size // 2)
    right_border = max_peak + (window_size // 2)
    if left_border < 0:
        shift = 0 - left_border
        left_border += shift
        right_border += shift
    elif right_border > max_x:
        shift = right_border - max_x
        left_border -= shift
        right_border -= shift

    return left_border, right_border


def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def generate_out_data(out_dir, directory_path, windows, file_2_pep_pos_dict, peptide_overlap):
    out_data = list()
    clust_2_file = {fasta_file.split('_')[-2]: fasta_file for fasta_file in sorted(windows.keys())}
    for cluster_id, fasta_file in clust_2_file.items():
        file_path = os.path.join(directory_path, fasta_file)

        fasta_dict = ft.read_fasta_dict(file_path)

        for pep_num, window in enumerate(windows[fasta_file], 1):
            for seq_name in sorted(fasta_dict.keys()):
                # check if original peptide is overlapping >= 9
                for probe_name, og_pep_window in file_2_pep_pos_dict[fasta_file].items():
                    og_seq_name = '_'.join(probe_name.split('_')[0:-2])
                    if og_seq_name == seq_name:
                        if get_overlap(og_pep_window, window) >= peptide_overlap:
                            new_pep_seq = fasta_dict[seq_name][window[0]:window[1]]
                            og_pep_seq = fasta_dict[seq_name][og_pep_window[0]:og_pep_window[1]]
                            out_data.append( (cluster_id, seq_name, f"Peptide_{pep_num}", new_pep_seq, probe_name, og_pep_seq) )

    out_df = pd.DataFrame(out_data, columns=["ClusterID", "SequenceName", "PeptideID", "NewPeptideSeq", "OriginalProbeName", "OriginalPeptideSeq"])
    out_df.to_csv(os.path.join(out_dir, "peptide_seq_data.tsv"), sep='\t', index=False)


def create_line_charts(alignCountsD, windows, out_dir):
    out_dir = os.path.join(out_dir, "final_window_visualizations")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for file, pos_dict in alignCountsD.items():
        create_line_chart(
                x = list(pos_dict.keys()),
                y = list(pos_dict.values()),
                windows = windows[file],
                title = file,
                out_dir = os.path.join(out_dir, f"{file.split('_')[-2]}_epitopes_lineplot.png")
                )

def create_line_chart(x, y, windows, out_dir, title):
    fig, ax = plt.subplots(figsize=(max(x)/10, 10), facecolor='w')
    ax.plot(x, y, linestyle='-')

    for window in windows:
        plt.axvspan(window[0], window[1], color="#ff6b0f", alpha=0.75)

    ax.set_xticks(np.arange(min(x), max(x)+5, 5))
    ax.set_xlim(left=min(x))
    ax.set_ylim(bottom=min(y))
    plt.grid()
    plt.xlabel("Sequence Position")
    plt.ylabel("Count")
    plt.title(title) 
    plt.savefig(out_dir, dpi=300, bbox_inches='tight')


def find_smallest_value_with_substring(data_dict, substring):
    # Filter the dictionary to only include items with the specified substring in the key
    filtered_dict = {k: v for k, v in data_dict.items() if substring in k}
    
    # If there are no matches, return None
    if not filtered_dict:
        return None
    
    # Find the key-value pair with the smallest value
    smallest_pair = min(filtered_dict.items(), key=lambda item: item[1])
    
    return smallest_pair


def read_check_align_file(directory):
    data_dict = {}
    clusters = set()

    # Construct the full file path
    filepath = os.path.join(directory, 'checkAlignLength.out')
    # Read the file content
    with open(filepath, 'r') as file:
        alignedCluster = None
        for line in file:
            if "mafft" in line:
                alignedCluster = line.strip()
                clusters.add(line.split('_')[-2])
            elif "Alignment:" in line and alignedCluster:
                alignLength = line.replace('Alignment:','').strip()
                #print(alignedCluster,alignLength)
                data_dict[alignedCluster] = alignLength
    # Find alignment with shortest length for each cluster
    results = {}
    for cluster in clusters:
        result = find_smallest_value_with_substring(data_dict, cluster)
        results[result[0]] = result[1]
                
    return results


def process_files_probes(probes_dict, directory_path):
    result = {}
    file_2_pep_pos_dict = dict()

    for filename, data in probes_dict.items():
        aligned_probes_file = filename.replace('.fasta', '_probesAligned.txt')
        aligned_probes_path = os.path.join(directory_path, aligned_probes_file)
        
        result[filename], file_2_pep_pos_dict[filename] = process_file_probes(data, aligned_probes_path)
    
    return result, file_2_pep_pos_dict

def process_file_probes(data, aligned_probes_path, removed_peptides = set()):
    aligned_length = int(data)
    
    alignD = {key: 0 for key in range(aligned_length + 1)}
    pep_pos_dict = dict()

    with open(aligned_probes_path, 'r') as file:
        for line_count, line in enumerate(file):
            if line_count > 0:
                elems = line.split('\t')
                if elems[0] not in removed_peptides:
                    seq_positions = elems[-1].split('~')
                    for pos in seq_positions:
                        alignD[int(pos)] += 1

                    pep_pos_dict[elems[0]] = (int(elems[1]), int(elems[2]))

    return alignD, pep_pos_dict


if __name__ == "__main__":
    main()
