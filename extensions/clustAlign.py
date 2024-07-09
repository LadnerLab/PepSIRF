#!/usr/bin/env python3
#-------DELETE LATER--------
import sys
sys.path.append("/Users/scg283/Documents/GitHub/Modules")
#---------------------------

import os
import glob
import pandas as pd
import fastatools as ft
import argparse
from collections import defaultdict
import subprocess
import concurrent.futures
import multiprocessing
import time

def main():
    start_time = time.perf_counter()

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i', '--bulk-clust-dir',  help='Directory created from bulk clustering script.', required=True)
    parser.add_argument('-q', '--quiet', required=False, action="store_true", help="Include to turn off mafft output messages.")
    parser.add_argument('--max-workers', type=int, required=False, default=None, help="Number of process to execute at a time.")
    parser.add_argument('-o', '--output-dir', default="out_dir", help='Name of Directory to output aligned and unaligned fasta files.')
    
    args = parser.parse_args()

    assert not os.path.exists(args.output_dir), \
        f"Warning: The directory \"{args.output_dir}\" already exists. Please move or delete."
    if args.max_workers != None:
        assert args.max_workers <= multiprocessing.cpu_count(), \
            f"Max workers excedes {multiprocessing.cpu_count()}, the number of CPUs on your machine."

    os.mkdir(args.output_dir)

    unaligned_dir =make_dir(args.output_dir, "unaligned")

    create_unaligned_dir(args.bulk_clust_dir, unaligned_dir)

    aligned_dir = make_dir(args.output_dir, "aligned")

    create_aligned_dir(unaligned_dir, aligned_dir, args.quiet, args.max_workers)

    align_pos_dir = make_dir(args.output_dir, "align_pos")

    create_align_pos_dir( aligned_dir, align_pos_dir )

    end_time = time.perf_counter()

    print(f"\nFinished in {round(end_time-start_time, 2)} seconds")


def create_align_pos_dir( aligned_dir, align_pos_dir ):
    for spec_dir in glob.glob(os.path.join(aligned_dir, '*')):
        out_species_dir = make_dir(align_pos_dir, os.path.basename(spec_dir))

        # loop through each thresh/multiprotein cluster directory
        for thresh_dir in glob.glob(os.path.join(spec_dir, '*')):
            out_thresh_dir = make_dir(out_species_dir, os.path.basename(thresh_dir))

            # loop through each protein/network directory
            for prot_dir in glob.glob(os.path.join(thresh_dir, '*')):
                out_prot_dir = make_dir( out_thresh_dir, os.path.basename(prot_dir))

                for in_fasta in glob.glob(os.path.join(prot_dir, '*')):
                    out_file = os.path.join(out_prot_dir, f"{in_fasta[len(os.path.join(prot_dir, '')):-len('.fasta')]}_positions.tsv")
                    create_align_pos_file(in_fasta, out_file)


def create_align_pos_file( in_fasta, out_file ):
    out_data = list()
    with open(in_fasta, 'r') as in_file:
        lines = in_file.read().splitlines()
        for line_idx in range(0, len(lines), 2):
            probe_name = lines[line_idx][1:]
            seq = lines[line_idx + 1]

            # get alignment positions of sequence
            align_pos = ""
            for pos, let in enumerate(seq, 1):
                if let != '-':
                    align_pos += f"{pos}~"
            align_pos = align_pos[:-1]

            out_data.append( (probe_name, align_pos) )

    pd.DataFrame(out_data, columns=["ProbeName", "AlignPos"]).to_csv(out_file, sep='\t', index=False)



def create_aligned_dir( unaligned_dir, aligned_dir, quiet, max_workers ):
    num_species = len(glob.glob(os.path.join(unaligned_dir, '*')))
    species_count = 0
    # process each species directory
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
        futures = [process_executor.submit( create_aligned_species_dir,
                                        aligned_dir, 
                                        spec_dir,
                                        quiet ) for spec_dir in glob.glob(os.path.join(unaligned_dir, '*'))]

        for future in concurrent.futures.as_completed(futures):
            species = future.result()
            species_count += 1
            print(f"{f'{species} complete.': <30} ({species_count}/{num_species})")


def create_aligned_species_dir( aligned_dir, spec_dir, quiet ):
    out_species_dir = make_dir(aligned_dir, os.path.basename(spec_dir))

    # loop through each thresh/multiprotein cluster directory
    for thresh_dir in glob.glob(os.path.join(spec_dir, '*')):
        out_thresh_dir = make_dir(out_species_dir, os.path.basename(thresh_dir))

        # loop through each protein/network directory
        for prot_dir in glob.glob(os.path.join(thresh_dir, '*')):
            out_prot_dir = make_dir( out_thresh_dir, os.path.basename(prot_dir))

            with concurrent.futures.ThreadPoolExecutor() as thread_executor:
                future = [thread_executor.submit( align_file,
                                                in_fasta, 
                                                os.path.join(out_prot_dir, os.path.basename(in_fasta)), 
                                                quiet) for in_fasta in glob.glob(os.path.join(prot_dir, '*'))]
            '''
            for in_fasta in glob.glob(os.path.join(prot_dir, '*')):
                align_file(in_fasta, os.path.join(out_prot_dir, os.path.basename(in_fasta)), quiet)
            '''
    return os.path.basename(spec_dir)


def align_file(in_fasta, out_fasta, quiet):
    command = ["mafft-einsi", in_fasta]
    if quiet:
        command.insert(1, "--quiet")

    result = subprocess.run(command, stdout = subprocess.PIPE, universal_newlines = True)

    with open(out_fasta, 'w') as out_fasta:
        # remove newline characters that mafft adds
        lines = result.stdout.split('\n')
        for line_idx in range(len(lines)):
            line = lines[line_idx]
            if line_idx == 0:
                out_fasta.write(line + '\n')
            elif line.startswith('>'):
                out_fasta.write('\n' + line + '\n')
            else:
                out_fasta.write(line)


def create_unaligned_dir( input_dir, output_dir ):
    # loop through each species
    for spec_dir in glob.glob(os.path.join(input_dir, '*')):
        if os.path.isdir(spec_dir) and len(glob.glob(os.path.join(spec_dir, "*.fasta"))) > 0:
            species = os.path.basename(spec_dir)
            out_species_dir = make_dir(output_dir, species)

            fasta_data = defaultdict(str)

            # create fasta dict
            for fasta_file in glob.glob(os.path.join(spec_dir, "*.fasta")):
                # merge fasta into dict
                fasta_data = {**fasta_data, **ft.read_fasta_dict(fasta_file)}

            # for each hierarchical cluster
            for clust_file in glob.glob(os.path.join(spec_dir, "hierarchical_clusters", "*.fasta.tsv")):
                # create directory for file
                prot = clust_file[len(os.path.join(spec_dir, "hierarchical_clusters", f"clusters_{species}_")):-len(".fasta.tsv")]

                # read file
                clustDf = pd.read_csv(clust_file, sep="\t").set_index("Sequence")

                # loop through each threshold
                for thresh in clustDf.columns.values.tolist():
                    create_unaligned_fasta_files(clustDf, thresh, fasta_data, out_species_dir, prot, species)

            # for each multiprotein network
            for net_dir in glob.glob(os.path.join(spec_dir, "multiprotein_networks", "network_*")):
                network = net_dir[len(os.path.join(spec_dir, "multiprotein_networks", '')):]

                #r read file
                clustDf = pd.read_csv(os.path.join(net_dir, f"{network}_multiprotein_cluster_sequences.tsv"), sep="\t").set_index("Sequence")

                create_unaligned_fasta_files(clustDf, clustDf.columns.values.tolist()[0], fasta_data, out_species_dir, network, species)


def create_unaligned_fasta_files(clustDf, col, fasta_data, out_species_dir, name, species):
    # make directories if they do not already exist
    col_out_dir = make_dir(out_species_dir, col)
    output_dir = make_dir(col_out_dir, f"{species}_{name}")

    clust_dict = create_clust_dict(clustDf, col)

    # create fasta files for each clust
    for clust, seqs in clust_dict.items():
        if len(seqs) > 1:
            # save to fasta file
            pd.DataFrame(
                [(f">{x}", fasta_data[x]) for x in seqs]
            ).to_csv(
                os.path.join(output_dir, f"{species}_{name}_{clust}.fasta"), sep="\n", index=False, header=False
            )


def make_dir( path, new ):
    dir_name = os.path.join(path, new)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    return dir_name


def create_clust_dict(df, thresh):
    clust_2_seq = defaultdict(list)
    seq_2_clust = df.to_dict()[thresh]

    for seq, clust in seq_2_clust.items():
        clust_2_seq[clust].append(seq)

    return clust_2_seq


if __name__ == "__main__":
    main()