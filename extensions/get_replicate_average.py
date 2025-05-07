#!/usr/bin/env python
import argparse


def create_list(current_line, item_to_remove):
    return_list = current_line.split("\t")
    if item_to_remove != "":
        return_list.remove(item_to_remove)
    return_list[len(return_list) - 1] = remove_newline(return_list[len(return_list) - 1])

    return return_list


def create_replicate_matrix(names_file, score_file, output_file):
    replicate_total = 0
    sequence_dict = {}
    base_sequences = []

    # Read in the name, score, and output files; print any errors found
    try:
        names_fh = open(names_file, "r")
        score_fh = open(score_file, "r")
        output_fh = open(output_file, "w")
    except FileNotFoundError as f:
        print(f"Error: file not found ({f})")
        return 1
    except NameError as n:
        print(f"Error: file not found ({n})")
        return 1
    except TypeError as t:
        print(f"Error: file not found ({t})")
        return 1

    output_fh.write("Sequence name\t")

    # Get the first line's base sequence
    current_line = names_fh.readline()
    current_line_list = create_list(current_line, "")
    base_sequence = current_line_list[0]

    # Write the first base sequence
    base_sequences.append(base_sequence)
    output_fh.write(f"{base_sequence}")
    sequence_dict[base_sequence] = []

    # Get the rest of the base sequences (if there are any)
    while True:
        # Get the current line
        current_line = names_fh.readline()
        if current_line == "":
            break

        # Get the base sequence on the current line
        current_line_list = create_list(current_line, "")
        base_sequence = current_line_list[0]

        # Use base sequence for the following:
        # 1. Append it to the list of base sequences
        # 2. Write it to the output file
        # 3. Add it as a key in the sequence dictionary
        base_sequences.append(base_sequence)
        output_fh.write("\t")
        output_fh.write(f"{base_sequence}")
        sequence_dict[base_sequence] = []

    output_fh.write("\n")

    # Get list of sequence names from input score file
    sequence_names_list = create_list(score_fh.readline(), "Sequence name")

    # Loop through score file & record the replicate averages
    while True:
        base_sequence_found = False

        # Check if at EOF
        current_line = score_fh.readline()
        if current_line == "":
            break

        # Create list out of values on current row
        current_line_list = create_list(current_line, "")

        # Get current peptide and remove from list
        current_peptide = current_line_list[0]
        current_line_list.pop(0)

        # Loop through sequence names
        sequence_names_index = 0
        while sequence_names_index < len(sequence_names_list):
            # Loop through base sequences
            base_sequence_index = 0
            while base_sequence_index < len(base_sequences):
                # Check if sequence name contains base sequence name
                if sequence_names_list[sequence_names_index].find(base_sequences[base_sequence_index]) != -1:
                    base_sequence_found = True

                    # Create/update sequence dictionary entry with:
                    # Key: sequence name
                    # Entry: list of replicate values
                    sequence_dict[base_sequences[base_sequence_index]].append\
                        (current_line_list[sequence_names_index])

                base_sequence_index += 1

            # Check if the base sequence was not found after the search
            if not base_sequence_found:
                print("Error: base sequence not found")

                names_fh.close()
                score_fh.close()
                output_fh.close()

                return 1

            sequence_names_index += 1

        output_fh.write(f"{current_peptide}")

        # For a given list of values at a sequence's entry,
        # Total up all the values and find their average
        for sequence in sequence_dict:
            for replicate in sequence_dict[sequence]:
                replicate_total += int(replicate)

            # Get the average and, if the average is an integer, write as is.
            # Otherwise, include precision for decimal value
            replicate_average = replicate_total / len(sequence_dict[sequence])
            output_fh.write("\t")
            if replicate_average.is_integer():
                output_fh.write(f"{int(replicate_average)}")
            else:
                output_fh.write('%.2f' % replicate_average)
            replicate_total = 0

        output_fh.write("\n")

        # Reset sequence dictionary
        for sequence in base_sequences:
            sequence_dict[sequence] = []

    # Close files
    names_fh.close()
    score_fh.close()
    output_fh.close()

    return 0


def main():
    arg_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg_parser.add_argument('-i', '--input', help="Name of input file that contains matrix with replicate totals")
    arg_parser.add_argument('-n', '--rep_names', help="Name of file that contains sequences to be utilized")
    arg_parser.add_argument('-a', '--get_avgs', help="Name of output file that matrix with replicate averages will go")

    args = arg_parser.parse_args()

    # Check for missing files; ensures that operations in create_replicate_matrix run correctly
    if args.rep_names is None and args.get_avgs is not None:
        print("Error: --get_avgs requires a file in --rep_names to be used")
        return 1
    elif args.get_avgs is None and args.rep_names is not None:
        print("Error: --rep_names requires a file in --get_avgs to be used")
        return 1

    create_replicate_matrix(args.rep_names, args.input, args.get_avgs)
    return 0


def remove_newline(string):
    edited_string = string
    if "\n" in string:
        edited_string = string.strip("\n")

    return edited_string


if __name__ == '__main__':
    main()
