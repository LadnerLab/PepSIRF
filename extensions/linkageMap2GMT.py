# sample command: python linkageMap2GMT.py -i PM1_linkageMap-species-k7_2024-03-26.tsv -o output.gmt
import argparse
import pandas as pd
import csv
import sys

parser = argparse.ArgumentParser(description="Convert Linkage Map into a GMT file")
parser.add_argument("-i", "--input-dir", type=str, metavar="", required=True, help="Filepath to input linkage map. [REQUIRED]")
parser.add_argument("-m", "--min-score", type=int, metavar="", required=False, default=1, help="Minimum score to filter species by. [DEFAULT=1]")
parser.add_argument("-o", "--output-dir", type=str, metavar="", required=True, help="Filepath where to write outputted GMT file. [REQUIRED]")
args=parser.parse_args()


# convert tsv linkage map to gmt
def linkage_map_2_gmt(
	input_filename: str,
	output_filename: str,
	min_score: int = 1
	)->None:
	
	df = pd.read_csv( 
		input_filename, sep="\t", header=0 
		).rename(columns={"Linked Species IDs with counts":"Species"}
		).dropna()

	# create dictionary with each peptide and its set of species
	peptide_dict = dict(zip(df["Peptide Name"], df["Species"]))

	for peptide, species in peptide_dict.items():
		peptide_dict[ peptide ] = str_to_set( str(species), min_score )

	# initialize create species dict
	species_dict = dict()
	for peptide in peptide_dict.keys():
		for species in peptide_dict[peptide]:
			if species not in species_dict.keys():
				species_dict[species] = peptide
			else:
				species_dict[species]+= f"\t{peptide}"

	with open( output_filename, "w" ) as gmt:
		for species in species_dict.keys():
			gmt.write(f"{species}\t\t{species_dict[species]}\n")

	print( f"Converted GMT file saved to: {output_filename}")

def str_to_set(
	val: str,
	min_score
	)->set:
	out_species = set();

	# get list of groups of species with counts
	spec_ls = val.split(",")

	for values in spec_ls:
		# species at 0 and score at 1
		split_id = values.split(":")

		# test if score greater or equal to than min
		if int(split_id[1]) >= min_score:
			# add all species to set (leave as string)
			out_species.add( split_id[0] )
			
	return out_species


if __name__ == "__main__":
	linkage_map_2_gmt( args.input_dir, args.output_dir, args.min_score )