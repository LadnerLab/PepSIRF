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
    parser.add_argument('-m', '--metadata',  help='Metadata file with info mapping the start and stop indexex of a sequence to a peptide name.', required=True)
    parser.add_argument('-q', '--quiet', required=False, action="store_true", help="Include to turn off mafft output messages.")
    parser.add_argument('--max-workers', type=int, required=False, default=None, help="Number of process to execute at a time.")
    parser.add_argument('-o', '--output-dir', default="out_dir", help='Name of Directory to output aligned and unaligned fasta files.')
    
    args = parser.parse_args()

    assert not os.path.exists(args.output_dir), \
        f"Warning: The directory \"{args.output_dir}\" already exists. Please move or delete."
    if args.max_workers != None:
        assert args.max_workers <= multiprocessing.cpu_count(), \
            f"Max workers excedes {multiprocessing.cpu_count()}, the number of CPUs on your machine."

    meta_df = pd.read_csv(args.metadata, sep='\t', usecols=["CodeName", "ParentSeqName", "StartIndex", "StopIndex"])

    os.mkdir(args.output_dir)

    unaligned_dir =make_dir(args.output_dir, "unaligned")

    create_unaligned_dir(args.bulk_clust_dir, unaligned_dir)

    aligned_dir = make_dir(args.output_dir, "aligned")

    create_aligned_dir(unaligned_dir, aligned_dir, args.quiet, args.max_workers, set(meta_df["ParentSeqName"]))

    align_pos_dir = make_dir(args.output_dir, "align_pos")

    create_align_pos_dir( aligned_dir, align_pos_dir, meta_df, args.max_workers )

    end_time = time.perf_counter()

    print(f"\nFinished in {round(end_time-start_time, 2)} seconds")


def create_align_pos_dir( aligned_dir, align_pos_dir, meta_df, max_workers ):
    for spec_dir in glob.glob(os.path.join(aligned_dir, '*')):
        out_species_dir = make_dir(align_pos_dir, os.path.basename(spec_dir))

        # loop through each network directory
        for net_dir in glob.glob(os.path.join(spec_dir, '*')):
            out_net_dir = make_dir(out_species_dir, os.path.basename(net_dir))

            # loop through each multiprotein cluster directory
            for prot_dir in glob.glob(os.path.join(net_dir, '*')):
                out_prot_dir = make_dir(out_net_dir, os.path.basename(prot_dir))

                # loop through each protein/network directory
                for mpc_dir in glob.glob(os.path.join(prot_dir, '*')):
                    out_mpc_dir = make_dir( out_prot_dir, os.path.basename(mpc_dir))

                    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        futures = [process_executor.submit( create_align_pos_file,
                                                        in_fasta, 
                                                        os.path.join(out_mpc_dir, f"{in_fasta[len(os.path.join(mpc_dir, '')):-len('.fasta')]}_positions.tsv"),
                                                        meta_df ) for in_fasta in glob.glob(os.path.join(mpc_dir, '*'))]


def create_align_pos_file( in_fasta, out_file, meta_df ):
    out_data = list()
    with open(in_fasta, 'r') as in_file:
        lines = in_file.read().splitlines()
        for line_idx in range(0, len(lines), 2):
            seq_name = lines[line_idx][1:]
            aligned_seq = lines[line_idx + 1]

            seq_df = meta_df.loc[meta_df["ParentSeqName"] == seq_name]

            pep_idx_dict = dict()
            for i, row in seq_df.iterrows():
                pep_idx_dict[row["CodeName"]] = (row["StartIndex"], row["StopIndex"])

            for pep, indices in pep_idx_dict.items():
                alignment = str()
                start = indices[0]
                stop = indices[1]

                seq_idx = 0
                aligned_seq_idx = 0

                while seq_idx < stop:
                    # safety for undefined behavior
                    if aligned_seq_idx >= len(aligned_seq):
                        raise IndexError("AAs in {seq_name} do not reach stop index for {pep}")
                    if aligned_seq[aligned_seq_idx] != '-':
                        if seq_idx >= start:
                            alignment += f"{aligned_seq_idx}~"
                        seq_idx += 1
                    aligned_seq_idx += 1

                out_data.append((pep, alignment[:-1]))

    pd.DataFrame(out_data, columns=["ProbeName", "AlignPos"]).to_csv(out_file, sep='\t', index=False)


def create_aligned_dir( unaligned_dir, aligned_dir, quiet, max_workers, seq_ids ):
    num_species = len(glob.glob(os.path.join(unaligned_dir, '*')))
    species_count = 0
    # process each species directory
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
        futures = [process_executor.submit( create_aligned_species_dir,
                                        aligned_dir, 
                                        spec_dir,
                                        quiet,
                                        seq_ids ) for spec_dir in glob.glob(os.path.join(unaligned_dir, '*'))]

        for future in concurrent.futures.as_completed(futures):
            species = future.result()
            species_count += 1
            print(f"{f'{species} alignment complete.': <30} ({species_count}/{num_species})")


def create_aligned_species_dir( aligned_dir, spec_dir, quiet, seq_ids ):
    out_species_dir = make_dir(aligned_dir, os.path.basename(spec_dir))

    # loop through each network directory
    for net_dir in glob.glob(os.path.join(spec_dir, '*')):
        out_net_dir = make_dir(out_species_dir, os.path.basename(net_dir))

        for prot_dir in glob.glob(os.path.join(net_dir, '*')):
            out_prot_dir = make_dir( out_net_dir, os.path.basename(prot_dir))

            for mpc_dir in glob.glob(os.path.join(prot_dir, '*')):
                out_mpc_dir = make_dir( out_prot_dir, os.path.basename(mpc_dir))

                for in_fasta in glob.glob(os.path.join(mpc_dir, '*')):
                    align_file(in_fasta, os.path.join(out_mpc_dir, os.path.basename(in_fasta)), quiet, seq_ids)

    return os.path.basename(spec_dir)


def align_file(in_fasta, out_fasta, quiet, seq_ids):
    command = ["mafft-einsi", "--preservecase", "--inputorder", "--thread", "-1", in_fasta]
    if quiet:
        command.insert(1, "--quiet")

    result = subprocess.run(command, stdout = subprocess.PIPE, universal_newlines = True)

    with open(out_fasta, 'w') as out_fasta:
        # remove newline characters that mafft adds
        lines = result.stdout.split('\n')
        for line_idx in range(len(lines)):
            line = lines[line_idx]

            # replace seq id if it is cut off by mafft
            if line.startswith('>'):
                if line[1:] not in seq_ids:
                    for seq_id in seq_ids:
                        if line[1:201] == seq_id[0:200]:
                            line = ">" + seq_id

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

            fasta_data = defaultdict(lambda: defaultdict(str))

            # create fasta dict
            for fasta_file in glob.glob(os.path.join(spec_dir, "*.fasta")):
                # merge fasta into dict
                fasta_data[fasta_file.split('_')[-1][:-len(".fasta")]] = ft.read_fasta_dict(fasta_file)

           
            # for each hierarchical cluster
            prot_2_df = defaultdict( lambda: pd.DataFrame )
            net_2_thresh = dict()
            for clust_file in glob.glob(os.path.join(spec_dir, "hierarchical_clusters", "*.fasta.tsv")):
                # create directory for file
                prot = clust_file[len(os.path.join(spec_dir, "hierarchical_clusters", f"clusters_{species}_")):-len(".fasta.tsv")]

                # read file
                prot_2_df[prot] = pd.read_csv(clust_file, sep="\t").set_index("Sequence")

                # loop through each threshold
                for net_num, thresh in enumerate(prot_2_df[prot].columns.values.tolist(), 1):
                    #create_unaligned_fasta_files(clustDf, thresh, fasta_data, out_species_dir, prot, species)
                    net_2_thresh[f"network_{net_num}"] = thresh
            

            # for each multiprotein network
            for net_dir in glob.glob(os.path.join(spec_dir, "multiprotein_networks", "network_*")):
                network = net_dir[len(os.path.join(spec_dir, "multiprotein_networks", '')):]

                # read file
                clustDf = pd.read_csv(os.path.join(net_dir, f"{network}_multiprotein_cluster_sequences.tsv"), sep="\t").set_index("Sequence")

                # make directories if they do not already exist
                col_out_dir = make_dir(out_species_dir, f"{network}")

                multi_clust_dict = create_clust_dict(clustDf, clustDf.columns.values.tolist()[0])

                # create fasta files for each clust
                for multi_clust, multi_clust_seqs in multi_clust_dict.items():
                    for prot in fasta_data:
                        # get all sequences from protein in this cluster 
                        prot_seqs = set(multi_clust_seqs).intersection(set(fasta_data[prot].keys()))
                        # check if there are sequences for that protein
                        if len(prot_seqs) > 0:
                            prot_dir = make_dir(col_out_dir, prot)

                            single_clust_dict = create_clust_dict_for_seqs(prot_2_df[prot], net_2_thresh[network], prot_seqs)

                            full_dir = make_dir(prot_dir, "full_mpc")

                            # check if multiple single protein clusters
                            if len(single_clust_dict.keys()) > 1:
                                indiv_dir = make_dir(prot_dir, "individual")
                                # output to individuals
                                for single_clust, single_clust_seqs in single_clust_dict.items():
                                    pd.DataFrame(
                                        [(f">{x}", fasta_data[prot][x]) for x in single_clust_seqs]
                                    ).to_csv(
                                        os.path.join(indiv_dir, f"mpc{multi_clust}_{single_clust}.fasta"), sep="\n", index=False, header=False
                                    )

                                # output concatenated to full mpc
                                pd.DataFrame(
                                    [(f">{seq}", fasta_data[prot][seq]) for seqs in single_clust_dict.values() for seq in seqs]
                                ).to_csv(
                                    os.path.join(full_dir, f"mpc{multi_clust}.fasta"), sep="\n", index=False, header=False
                                )
                            else:
                                # output concatenated to full mpc
                                pd.DataFrame(
                                    [(f">{seq}", fasta_data[prot][seq]) for seqs in single_clust_dict.values() for seq in seqs]
                                ).to_csv(
                                    os.path.join(full_dir, f"mpc{multi_clust}.fasta"), sep="\n", index=False, header=False
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

def create_clust_dict_for_seqs(df, thresh, seqs):
    clust_2_seq = defaultdict(list)
    seq_2_clust = {seq: clust for seq, clust in df.to_dict()[thresh].items() if seq in seqs}

    for seq, clust in seq_2_clust.items():
        clust_2_seq[clust].append(seq)

    return clust_2_seq


if __name__ == "__main__":
    main()