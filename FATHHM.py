# FATHHM.py
# Merges FATHMM pathogenicity scores with variant files based on position and reference base

import os
import pandas as pd

def FATHMM_merge(df1, file, i):
    """
    Merges FATHMM scores with a variant file.

    Args:
        df1 (pd.DataFrame): FATHMM data (tab-delimited, with position info).
        file (pd.DataFrame): Variant CSV file.
        i (int): Sample/group ID to be used in output filename.
    """
    df1.drop([0], inplace=True)  # Drop header row if duplicated
    df1.reset_index(drop=True, inplace=True)
    df1['Position'] = df1['Position'].astype(int)

    # Merge on position and reference base
    merged = pd.merge(file, df1, how='left', left_on=['left', 'ref_seq'], right_on=['Position', 'Ref. Base'])

    # Drop redundant columns after merge
    merged = merged.drop(columns=['# Chromosome', 'Position', 'Ref. Base', 'Mutant Base'])

    # Save merged result to file
    merged.to_csv(f"{i}_common_combine.csv", index=False)
    return

if __name__ == '__main__':
    # Set directory containing FATHMM and variant files
    os.chdir("./data/fathmm_input")

    
    fathmm_df = pd.read_csv("3_fathmm_common.txt", sep='\t')
    variant_df = pd.read_csv("3_common_combine.csv")

    FATHMM_merge(fathmm_df, variant_df, 3)
    print("FATHMM merge complete.")
