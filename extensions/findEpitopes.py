import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i', '--input-dir',  help='Directory with alignment output files and files that contain the mapped location of peptides', required=True)
    parser.add_argument('--window-size', type=int, default=30,  help='Size of AA window to use for identifying core epitopes.', required=False)
    parser.add_argument('--max-zeros', type=int, default=5, help='Maximum number of zero counts a window can contain.', required=False)
    parser.add_argument('--max-overlap', type=int, default=8, help='Maximum AA overlap a window can have with a previously selected window.', required=False)
    parser.add_argument('-o', '--output-dir', default="clust_align_visualizations", help='Name of directory to output line plots.')
    
    args = parser.parse_args()

    directory_path = args.input_dir
    alignment_to_use_dict = read_check_align_file(directory_path)
    #print(probes_dict)
    alignCountsD = process_probes(alignment_to_use_dict, directory_path)
    #print(alignCountsD)

    windows = find_core_epitopes(alignCountsD, args.window_size, args.max_zeros, args.max_overlap)
    print(windows)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    create_line_chart(alignCountsD, windows, args.output_dir)


def find_core_epitopes(alignCountsD, window_size, max_zeros, max_overlap):
    out_dict = dict()

    for file in alignCountsD.keys():
        windows = dict()

        counts = list(alignCountsD[file].values())

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

                # check no more than max zeros and does not overlap any previously selected window by more than max overlap
                if window.count(0) <= max_zeros and all(get_overlap((start_idx, end_idx), (x[0], x[1])) <= max_overlap for x in list(windows.keys())):
                    # check if greater than max score
                    score = sum(window)                    
                    if score > max_score:

                        # center window around peak
                        possible_windows = [(start_idx, end_idx)]
                        temp_start = start_idx + 1
                        temp_end = end_idx + 1
                        temp_window = counts[temp_start:temp_end]
                        while sum(temp_window) == score and temp_window.count(0) <= max_zeros and all(get_overlap((temp_start, temp_end), (x[0], x[1])) <= max_overlap for x in list(windows.keys())):
                            possible_windows.append((temp_start,temp_end))
                            temp_start += 1
                            temp_end += 1
                            temp_window = counts[temp_start:temp_end] 

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
                windows[max_window] = max_score

        out_dict[file] = windows

    return out_dict


def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def create_line_chart(alignCountsD, windows, out_dir):
    for file, pos_dict in alignCountsD.items():
        x = list(pos_dict.keys())
        y = list(pos_dict.values())

        fig, ax = plt.subplots(figsize=(max(x)/10, 10), facecolor='w')
        ax.plot(x, y, linestyle='-')
        #cmap = mpl.colormaps['Oranges']
        #norm = colors.Normalize(vmin=0, vmax=len(windows[file])-1)
        for idx, window in enumerate(list(windows[file].keys())):
            plt.axvspan(window[0], window[1], color="#ff6b0f", alpha=0.75)

        ax.set_xticks(np.arange(min(x), max(x)+5, 5))
        ax.set_xlim(left=min(x))
        ax.set_ylim(bottom=min(y))
        plt.grid()
        plt.xlabel("Sequence Position")
        plt.ylabel("Count")
        plt.title(file) 
        plt.savefig(os.path.join(out_dir, f"{file.split('_')[-2]}_epitopes_lineplot.png"), dpi=300, bbox_inches='tight')


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


def process_probes(probes_dict, directory_path):
    result = {}

    for filename, data in probes_dict.items():
        aligned_probes_file = filename.replace('.fasta', '_probesAligned.txt')
        aligned_probes_path = os.path.join(directory_path, aligned_probes_file)
        
        aligned_length = int(data)
        #print(range(0, aligned_length))
        
        alignD = {key: 0 for key in range(aligned_length + 1)}

        with open(aligned_probes_path, 'r') as file:
            for line_count, line in enumerate(file):
                if line_count > 0:
                    seq_positions = line.split('\t')[-1].split('~')
                    for pos in seq_positions:
                        alignD[int(pos)] += 1
        
        result[filename] = alignD
    
    return result

if __name__ == "__main__":
    main()
