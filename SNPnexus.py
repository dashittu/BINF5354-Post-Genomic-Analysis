# SNPnexus.py
# Prepares SNP input files for SNPnexus and merges SNPnexus results with variant files

import os
import pandas as pd

def SNPnexus_input(csvfile, n):
    """
    Extracts relevant information from a CSV file and creates a text input file
    for SNPnexus annotation.

    Args:
        csvfile (pd.DataFrame): Input variant data
        n (int): Identifier to use in output filename
    """
    pos_variant = []
    ref = []
    alt = []
    chrom = []
    strand = []
    num_chrom = []

    for line in range(csvfile.shape[0]):
        chrom.append("Chromosome")
        num_chrom.append(csvfile["chrom"][line].replace("chr", ""))
        pos_variant.append(csvfile["left"][line])
        ref.append(csvfile["ref_seq"][line])
        strand.append(1)

        if csvfile["ref_seq"][line] != csvfile["var_seq1"][line]:
            alt.append(csvfile["var_seq1"][line])
        else:
            alt.append(csvfile["var_seq2"][line])

    # Create DataFrame for SNPnexus format
    dp = pd.DataFrame(chrom, columns=['chromosome'])
    dp.insert(1, "num_chrom", num_chrom)
    dp.insert(2, "position", pos_variant)
    dp.insert(3, "ref", ref)
    dp.insert(4, "alt", alt)
    dp.insert(5, "strand", strand)

    # Output to tab-separated txt file
    dp.to_csv(f"{n}_SNPnexus_input_tumor.txt", sep="\t", index=False)
    return

def SNPnexus_merge(df1, df2, file, i):
    """
    Merges SNPnexus annotation results into the original variant CSV file.

    Args:
        df1 (pd.DataFrame): SNPnexus near gene output
        df2 (pd.DataFrame): SNPnexus genomic coordinate output
        file (pd.DataFrame): Original variant file to be annotated
        i (int): Identifier to use in output filename
    """
    file = pd.merge(file, df1, how='left', left_on=['chrom', 'left'], right_on=['Chromosome', 'Position'])
    file = file.drop(columns=['Variation ID', 'Chromosome', 'Position', 'Type', 'Nearest Upstream Gene',
                              'Type of Nearest Upstream Gene', 'Distance to Nearest Upstream Gene',
                              'Nearest Downstream Gene', 'Type of Nearest Downstream Gene',
                              'Distance to Nearest Downstream Gene'])

    file = pd.merge(file, df2, how='left', left_on=['left'], right_on=['Position'])
    file = file.drop(columns=['Variation ID', 'Chromosome', 'Position', 'REF Allele', 'ALT Allele (IUPAC)',
                              'Minor Allele', 'Minor Allele Global Frequency', 'Contig',
                              'Contig Position', 'Band'])

    file.to_csv(f"{i}_common_combine.csv", index=False)
    return

if __name__ == '__main__':
    # --------- Step 1: Generate SNPnexus input files ---------
    os.chdir("./data/SNPnexus_input")  # Update path to your local SNPnexus input folder
    csv_files = os.listdir()

    n = 1
    for file in csv_files:
        if file.endswith(".csv"):
            SNPnexus_input(pd.read_csv(file), n)
            n += 1

    # --------- Step 2: Merge SNPnexus annotations back ---------
    # Note: Update file names for your specific data
    df_near = pd.read_csv("3_near_gens.txt", sep='\t')
    df_coords = pd.read_csv("3_gen_coords.txt", sep='\t')
    base_file = pd.read_csv("3_common.csv")

    SNPnexus_merge(df_near, df_coords, base_file, 3)
    print("SNPnexus input and merge complete.")
