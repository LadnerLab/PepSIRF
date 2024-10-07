import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
import fastatools as ft

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

    generate_out_data(args.output_dir, directory_path, last_iter_counts, alignCountsD, file_2_pep_pos_dict, windows, args.window_size, args.peak_overlap_window_size)

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
        while peak_found:
            peak_found = False
            alignCountsD, pep_pos_dict = process_file_probes(
                                        data=data, 
                                        aligned_probes_path=aligned_probes_path,
                                        removed_peptides=removed_peptides
                                            )
            max_window, number_of_windows = find_core_epitopes(alignCountsD, window_size, max_zeros, max_overlap, windows)
            if number_of_windows > 0:
                peak_found = True
                left_border,right_border = max_window
                pep_ovlp_win = (max_window[0] + 10, max_window[1] - 10)
                
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


def generate_out_data(out_dir, directory_path, last_iter_counts, alignCountsD, pep_pos_dict, windows, window_size, peak_ovlp_window_size):
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
            pep_ovlp_win = (window[0] + 10, window[1] - 10)
            # Count overlapping peptides that completely middle 10AA cover window
            for pep, pep_coords in pep_pos_dict[fasta_file].items():
                if get_overlap(pep_ovlp_win, pep_coords) == peak_ovlp_window_size:
                    ovlp_count += 1

            out_data.append( (fasta_file, cluster_id, f"Peptide_{pep_num}", window[0], window[1], round(counts_avg, 3), counts_sum, ovlp_count) )

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

    out_df = pd.DataFrame(out_data, columns=["Fasta File","ClusterID", "PeptideID", "Start Position", "Stop Position", "Peptide Counts Average", "Peptide Counts Sum", "Peptide Overlap"])
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
    width = max(x)/10
    if width > 100:
        width = 100
    fig, ax = plt.subplots(figsize=(width, 10), facecolor='w')
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
    #filepath = os.path.join('checkAlignLength.out')

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
    counts = list(alignCountsD.values())
    possible_windowsD = {}
    possible_windows_avgD = {}
    start_idx = 0
    while start_idx < len(counts) - window_size + 1:
        end_idx = start_idx + window_size
        window = counts[start_idx:end_idx]

        # Check no more than max_zeros and no overlap beyond max_overlap
        if window.count(0) <= max_zeros and all(get_overlap((start_idx, end_idx), x) <= max_overlap for x in peak_windows):
            possible_windowsD[(start_idx, end_idx)] = sum(window)
            possible_windows_avgD[(start_idx, end_idx)] = sum(window) / len(window)

        start_idx += 1

    number_of_windows = len(possible_windowsD)
    if number_of_windows == 0:
        return None, 0
    # If remaining possible windows have an average score > 1.3333 (overlap of 10 or more
    # positions for single overlapping peptides)
    if max(possible_windows_avgD.values()) >= 1.3333:
        # Find max score window(s)
        max_score = max(possible_windowsD.values())
        max_windowList = [w for w, score in possible_windowsD.items() if score >= max_score * 0.9]
        print(f"max_windowList: {max_windowList}")
        for ele in max_windowList:
            print(ele,possible_windowsD[ele])
    
        # Find window with highest center score
        mid_window_scoresD = {}
        max_mid_window_score = 0
        for window in max_windowList:
            mid_window_score = sum(counts[window[0] + (window_size // 3):window[1] - (window_size // 3)])
            mid_window_scoresD[window] = mid_window_score
            if mid_window_score > max_mid_window_score:
                max_mid_window_score = mid_window_score
                max_window = window
        print(f"mid_window_scoresD: {mid_window_scoresD}")
        # If multiple windows have the same middle window score, pick the window with the largest total score
        max_mid_windows = [w for w, score in mid_window_scoresD.items() if score == max_mid_window_score]
        if len(max_mid_windows) > 1:
            max_window = max(max_mid_windows, key=lambda w: possible_windowsD[w])
            max_window_score = possible_windowsD[max_window]
            max_window_scores = [w for w, score in possible_windowsD.items() if score == max_window_score]
            max_full_mid_window_scores = [w for w in max_window_scores if w in max_mid_windows]
            # If there are multiple windows with the same full window score, take the middle one
            if len(max_full_mid_window_scores) > 1:
                max_window = max_full_mid_window_scores[len(max_full_mid_window_scores)//2]  
                print(max_window)
    # If remaining possible windows are all individual peptides with little overlap
    # switch to simple sliding window approach choosing all windows that pass thresholds
    if  max(possible_windows_avgD.values()) < 1.3333:
        start_count = alignCountsD[list(possible_windowsD.keys())[0][0]]
        possible_windows_noZeroStart = [w for w, score in possible_windowsD.items() if alignCountsD[w[0]] != 0]
        print(f"possible_windows_noZeroStart: {possible_windows_noZeroStart}")
        max_window = possible_windows_noZeroStart[0]
        print(max_window)
    
    return max_window, number_of_windows

if __name__ == "__main__":
    main()
