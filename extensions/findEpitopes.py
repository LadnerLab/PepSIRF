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

    windows, last_iter_counts = iterative_peptide_finder(alignment_to_use_dict, directory_path, args.window_size, args.max_zeros, args.max_overlap, 
                                                                args.peak_overlap_window_size, args.output_dir, args.include_iter_vis)

    # generate_out_data(args.output_dir, directory_path, windows, file_2_pep_pos_dict, args.peptide_overlap)

    generate_out_data(args.output_dir, directory_path, last_iter_counts, alignCountsD, windows, args.window_size)

    create_line_charts(alignCountsD, windows, args.output_dir)


def iterative_peptide_finder(alignment_to_use_dict, directory_path, window_size, max_zeros, max_overlap, 
                                                                    peak_ovlp_window_size, output_dir, include_iter_vis):
    out_dict = dict()
    last_iter_counts = dict()

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
        invalid_peaks_found = False
        while peak_found:
            if not invalid_peaks_found:
                peak_found = False
                score = 0

                # generate scores without used peptides
                alignCountsD, pep_pos_dict = process_file_probes(
                                                        data=data, 
                                                        aligned_probes_path=aligned_probes_path,
                                                        removed_peptides=removed_peptides
                                                        )

                y = list(alignCountsD.values())
                
                # get all peaks
                _, peak_plateaus = find_peaks(y, plateau_size = 1)
                peaks = [(peak_plateaus["left_edges"][i], peak_plateaus["right_edges"][i]) for i in range(len(peak_plateaus["plateau_sizes"]))]
                
                if len(peaks) > 0:
                    valid_window = False
                    valid_peaks_exist = True
                    invalid_peaks = set()
                    while not valid_window and valid_peaks_exist:
                        # get max borders (such that entire plateau is not invalid)
                        valid_peaks = [peak for peak in peaks if any(x not in invalid_peaks for x in range(peak[0], peak[1] + 1))]
                        if valid_peaks:
                            max_peak_borders_ties = list()
                            max_height = 0
                            for peak in valid_peaks:
                                if y[peak[0]] > max_height:
                                    max_peak_borders_ties.clear()
                                    max_height = y[peak[0]]
                                    max_peak_borders_ties.append(peak)
                                elif y[peak[0]] == max_height:
                                    max_peak_borders_ties.append(peak)


                            # find which point has the highest score across the window
                            max_peak_window_score = 0
                            max_peak_ties=list()
                            for max_peak_borders in max_peak_borders_ties:
                                for peak in range(max_peak_borders[0], max_peak_borders[1] + 1):
                                    if peak not in invalid_peaks:
                                        # center window on whichever has greater sum
                                        left_border_l, right_border_l = generate_window(peak, window_size, int(data), "left")
                                        left_border_r, right_border_r = generate_window(peak, window_size, int(data), "right")
                                        if sum(y[left_border_l:right_border_l]) >= sum(y[left_border_r:right_border_r]) or window_size % 2 == 0:
                                            left_border = left_border_l
                                            right_border = right_border_l
                                        else:
                                            left_border = left_border_r
                                            right_border = right_border_r

                                        # generate window will 
                                        total_window_score = sum(y[left_border:right_border])
                                        if total_window_score > max_peak_window_score:
                                            max_peak_ties.clear()
                                            max_peak_window_score = total_window_score
                                            max_peak_ties.append(peak)
                                        elif total_window_score == max_peak_window_score:
                                            max_peak_ties.append(peak)

                            max_peak = max_peak_ties[len(max_peak_ties)//2]

                            # get designed peptide window
                            left_border, right_border = generate_window(max_peak, window_size, int(data))

                            # test if window passes thresholds
                            if y[left_border:right_border].count(0) <= max_zeros and \
                                            all(get_overlap((left_border, right_border), (x[0], x[1])) <= max_overlap for x in windows):
                                valid_window = True
                            else:
                                # remove possible mass peak
                                invalid_peaks.add(max_peak)
                                invalid_peaks_found = True
                                    
                        else:
                            valid_peaks_exist = False

                    if valid_peaks_exist:
                        peak_found = True
                        pep_ovlp_win = generate_window(max_peak, peak_ovlp_window_size, int(data))

                        # remove overlapping peptides that completely cover window
                        for pep, pep_coords in pep_pos_dict.items():
                            if get_overlap(pep_ovlp_win, pep_coords) == peak_ovlp_window_size:
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
                    else:
                        peak_found = False

            else:
                alignCountsD, pep_pos_dict = process_file_probes(
                                            data=data, 
                                            aligned_probes_path=aligned_probes_path,
                                            removed_peptides=removed_peptides,
                                                )
                max_window, number_of_windows = find_core_epitopes(alignCountsD, window_size, max_zeros, max_overlap, windows)
                if number_of_windows == 0:
                    break
                if number_of_windows == 1:
                    peak_found = False
                    invalid_peaks_found = False
                left_border,right_border = max_window
                # remove overlapping peptides that completely cover window
                for pep, pep_coords in pep_pos_dict.items():
                    if get_overlap(max_window, pep_coords) == peak_ovlp_window_size:
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
        last_iter_counts[filename] = alignCountsD

    return out_dict, last_iter_counts


def generate_window(max_peak, window_size, max_x, how="left"):
    left_border = max_peak - (window_size // 2)
    right_border = max_peak + (window_size // 2)

    if how == "right":
        left_border += 1
        right_border += 1
    elif how != "left":
        raise NotImplementedError(f"{how} not a valid argumentd for how to generate window.")

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


def generate_out_data(out_dir, directory_path, last_iter_counts, alignCountsD, windows, window_size):
    out_data = list()
    sum_data = list()
    clust_2_file = {fasta_file.split('_')[-2]: fasta_file for fasta_file in sorted(windows.keys())}
    for cluster_id, fasta_file in clust_2_file.items():
        # create output file
        file_path = os.path.join(directory_path, fasta_file)

        fasta_dict = ft.read_fasta_dict(file_path)

        for pep_num, window in enumerate(windows[fasta_file], 1):
            counts_sum = sum(list(alignCountsD[fasta_file].values())[window[0]:window[1]])
            counts_avg = counts_sum / window_size

            ovlp_count = 0
            for pep_num_2, window_2 in enumerate(windows[fasta_file], 1):
                if pep_num != pep_num_2:
                    if get_overlap(window, window_2) > 0:
                        ovlp_count += 1

            out_data.append( (cluster_id, f"Peptide_{pep_num}", window[0], window[1], round(counts_avg, 3), counts_sum, ovlp_count) )

        # create summary stats
        not_covered_pos_count = 0
        not_covered_total = 0
        for pos, count in last_iter_counts[fasta_file].items():
            # check is position is not covered by any window
            if count > 0 and all([pos not in range(window[0],window[1]) for window in windows[fasta_file]]):
                not_covered_pos_count += 1
                not_covered_total += count

            norm_not_covered_pos_count = not_covered_pos_count / max([int(pos) for pos in last_iter_counts[fasta_file].keys()])

        sum_data.append( (cluster_id, not_covered_pos_count, round(norm_not_covered_pos_count, 3), not_covered_total ) )

    out_df = pd.DataFrame(out_data, columns=["ClusterID", "PeptideID", "Start Position", "Stop Position", "Peptide Counts Average", "Peptide Counts Sum", "Peptide Overlap"])
    out_df.to_csv(os.path.join(out_dir, "peptide_seq_data.tsv"), sep='\t', index=False)

    sum_df = pd.DataFrame(sum_data, columns=["ClusterID", "Number of Uncovered Positions (with count > 0)", "Normalized Number of Uncovered Positions", "Uncovered Positions Counts Sum"])
    sum_df.to_csv(os.path.join(out_dir, "summary_data.tsv"), sep='\t', index=False)


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
        plt.axvspan(window[0], window[1] - 1, color="#ff6b0f", alpha=0.75)

    ax.set_xticks(np.arange(min(x), max(x)+5, 5))
    ax.set_xlim(left=min(x))
    ax.set_ylim(bottom=min(y))
    plt.grid()
    plt.xlabel("Sequence Position")
    plt.ylabel("Count")
    plt.title(title) 
    plt.savefig(out_dir, dpi=300, bbox_inches='tight')
    plt.close()


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
    

def find_core_epitopes(alignCountsD, window_size, max_zeros, max_overlap, peak_windows):
    peak_windows = peak_windows
    new_windows = dict()

    counts = list(alignCountsD.values())

    # get windows
    window_found = True
    while window_found:
        window_found = False
        max_score = 0
        # iterate through each possible window
        start_idx = 0
        while start_idx < len(counts) - window_size + 1:
            end_idx = start_idx + window_size
            window = counts[start_idx:end_idx]
            current_windows = peak_windows + list(new_windows.keys())

            # check no more than max zeros and does not overlap any previously selected window by more than max overlap
            if window.count(0) <= max_zeros and all(get_overlap((start_idx, end_idx), (x[0], x[1])) <= max_overlap for x in current_windows):
                # check if greater than max score
                score = sum(window)                    
                if score > max_score:

                    # center window around peak
                    possible_windows = [(start_idx, end_idx)]
                    temp_start = start_idx + 1
                    temp_end = end_idx + 1
                    temp_window = counts[temp_start:temp_end]
                    while sum(temp_window) == score and temp_window.count(0) <= max_zeros and all(get_overlap((temp_start, temp_end), (x[0], x[1])) <= max_overlap for x in current_windows):
                        possible_windows.append((temp_start,temp_end))

                        # increment window
                        temp_start += 1
                        temp_end += 1
                        temp_window = counts[temp_start:temp_end] 

                    # find window with highest center score
                    mid_window_scores = dict()
                    for window in possible_windows:
                        mid_window_scores[window] = sum(counts[window[0] + (window_size // 3):window[1] - (window_size // 3)])
                    possible_windows = [x for x, y in mid_window_scores.items() if y == max(list(mid_window_scores.values()))]

                    # use smaller median from odd cases
                    if len(possible_windows) % 2 == 0:
                        start_idx, end_idx = possible_windows[(len(possible_windows) // 2) - 1]
                    else:
                        start_idx, end_idx = possible_windows[(len(possible_windows) // 2)]

                    # set max
                    max_score = score
                    max_window = (start_idx, end_idx)
                    window_found = True

            start_idx += 1

        if window_found:
            new_windows[max_window] = max_score
                
    number_of_windows = len(new_windows)
    
    if number_of_windows > 0:
        max_window = max(zip(new_windows.values(), new_windows.keys()))[1]
    if number_of_windows == 0:
        max_window = None

    return max_window, number_of_windows


if __name__ == "__main__":
    main()
