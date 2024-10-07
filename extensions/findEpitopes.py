import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i', '--input-dir',  help='Directory with alignment output files and files that contain the mapped location of peptides', required=True)
    parser.add_argument('-o', '--output-dir', default="clust_align_visualizations", help='Name of directory to output line plots.')
    
    args = parser.parse_args()

    directory_path = args.input_dir
    alignment_to_use_dict = read_check_align_file(directory_path)
    #print(probes_dict)
    alignCountsD = process_probes(alignment_to_use_dict, directory_path)
    #print(alignCountsD)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    create_line_chart(alignCountsD, args.output_dir)


def create_line_chart(alignCountsD, out_dir):
    for file, pos_dict in alignCountsD.items():
        x = list(pos_dict.keys())
        y = list(pos_dict.values())
        fig, ax = plt.subplots(figsize=(max(x)/10, 10), facecolor='w')
        ax.plot(x, y, linestyle='-')
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
