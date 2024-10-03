import os
import argparse
import pandas as pd
import math

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    parser.add_argument('-i', '--input-filepath', type=str, help='Path to tab-separated file with columns for sample, concentration, and volume', required=True)
    parser.add_argument('--max-vol-sum', type=float, default=None, help='Maximum sum of all of the volumes. Set to sum of volume avaiable by default.', required=False)
    parser.add_argument('-o', '--output-filepath', type=str, default="dilutions_out.tsv", help='Filepath to output final dilutions', required=False)
    
    args = parser.parse_args()

    input_df = pd.read_csv(args.input_filepath, sep='\t')

    final_dilutions_df = calculate_dilutions(
                input_df = input_df,
                max_vol_sum = args.max_vol_sum
    )

    final_dilutions_df.to_csv(args.output_filepath, sep='\t', index=False)


def calculate_dilutions(
                input_df, 
                max_vol_sum,
) -> pd.DataFrame:
    # set no min vol sum if not given
    current_vol_sum = sum(input_df["vol_available"])
    if max_vol_sum:
        assert max_vol_sum <= current_vol_sum, \
            f"--max_vol_sum must be less or equal than the current volume sum, {current_vol_sum}"
    else:
        max_vol_sum = current_vol_sum

    all_c1 = input_df["conc"].to_list()
    all_vol_available = input_df["vol_available"].to_list()

    # get c2
    c2 = 1 / ( sum( [ 1 / c1 for c1 in all_c1 ]))

    v2 = max_vol_sum
    all_v1 = list()
    valid_volumes = False
    while not valid_volumes:
        for idx in range(len(all_c1)):
            c1 = all_c1[ idx ]
            vol_available = all_vol_available[ idx ]

            v1 = (c2 / c1) * v2

            # check if v1 is invalid
            if v1 > vol_available:
                all_v1.clear()
                v2 -= 0.1
                break
            else:
                all_v1.append(v1)

        # check if for loop ended
        if len(all_v1) == len(all_c1):
            valid_volumes = True

    print(f"Volume sum: {round(v2, 1)}")

    final_dilutions_df = pd.DataFrame(([(round(c2, 3), round(v1, 2)) for v1 in all_v1]), columns=["conc_pool", "vol_pool"])
    final_dilutions_df = input_df.join(final_dilutions_df)
    return final_dilutions_df


if __name__ == "__main__":
    main()