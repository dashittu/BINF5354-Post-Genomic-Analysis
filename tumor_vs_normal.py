# tumor_vs_normal.py
# Compares tumor vs normal variant files, identifies shared and unique variants,
# and generates subject-level mutation statistics.

import os
import pandas as pd

# Set directories (update as needed)
tumordir = "./data/tumor/"
normaldir = "./data/normal/"

# Get file names
tumors = [f for f in os.listdir(tumordir) if f.endswith('.csv')]
normals = [f for f in os.listdir(normaldir) if f.endswith('.csv')]
tumorname = sorted([f.split("_")[0] for f in tumors])
normalname = sorted([f.split("_")[0] for f in normals])

# Only proceed if tumor and normal filenames match
if tumorname == normalname:
    disease_var = pd.DataFrame()
    normal_var = pd.DataFrame()
    common_var = pd.DataFrame()

    for name in tumorname:
        tf = pd.read_csv(tumordir + name + "_tumor.csv")
        nf = pd.read_csv(normaldir + name + "_normal.csv")

        # Drop duplicate variants within individual files
        norepeat_tf = tf.drop_duplicates(subset=['chrom', 'left', 'ref_seq', 'var_seq1', 'var_seq2'], keep='first')
        norepeat_nf = nf.drop_duplicates(subset=['chrom', 'left', 'ref_seq', 'var_seq1', 'var_seq2'], keep='first')

        # Accumulate variant records
        disease_var = disease_var.append(norepeat_tf, ignore_index=True)
        normal_var = normal_var.append(norepeat_nf, ignore_index=True)

    # Filter out rows with 'NoName' gene
    disease_var = disease_var[disease_var['genename'] != 'NoName'].reset_index(drop=True)
    normal_var = normal_var[normal_var['genename'] != 'NoName'].reset_index(drop=True)

    # Identify shared variants
    common_var = pd.merge(disease_var, normal_var, on=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2'], how='inner')
    common_var = common_var.drop(columns=['var_index_y', 'var_score_y', 'genename_y', 'where_y', 'change_type1_y'])
    common_var.rename(columns={
        'var_index_x': 'var_index',
        'var_score_x': 'var_score',
        'genename_x': 'genename',
        'where_x': 'where',
        'change_type1_x': 'change_type1'
    }, inplace=True)

    # Subtract shared variants from disease and normal to get unique sets
    disease_var = disease_var.append(common_var).drop_duplicates(keep=False, subset=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2'])
    normal_var = normal_var.append(common_var).drop_duplicates(keep=False, subset=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2'])

    # Drop unnecessary columns
    disease_var.drop(["var_index", "var_score"], axis=1, inplace=True)
    normal_var.drop(["var_index", "var_score"], axis=1, inplace=True)
    common_var.drop(["var_index", "var_score"], axis=1, inplace=True)

    # Compute subject-level mutation counts for each variant group
    disease_all = pd.DataFrame()
    normal_all = pd.DataFrame()
    common_all = pd.DataFrame()

    for name in tumorname:
        tf = pd.read_csv(tumordir + name + "_tumor.csv")
        nf = pd.read_csv(normaldir + name + "_normal.csv")

        norepeat_tf = tf.drop_duplicates(subset=['chrom', 'left', 'ref_seq', 'var_seq1', 'var_seq2'])
        norepeat_nf = nf.drop_duplicates(subset=['chrom', 'left', 'ref_seq', 'var_seq1', 'var_seq2'])

        common = pd.merge(norepeat_tf, norepeat_nf, on=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2'], how='inner')
        common = common.drop(columns=['var_index_y', 'var_score_y', 'genename_y', 'where_y', 'change_type1_y'])
        common.rename(columns={
            'var_index_x': 'var_index',
            'var_score_x': 'var_score',
            'genename_x': 'genename',
            'where_x': 'where',
            'change_type1_x': 'change_type1'
        }, inplace=True)

        if not common.empty:
            norepeat_tf = norepeat_tf.append(common).drop_duplicates(keep=False, subset=['chrom', 'left', 'ref_seq', 'var_seq1', 'var_seq2'])
            norepeat_nf = norepeat_nf.append(common).drop_duplicates(keep=False, subset=['chrom', 'left', 'ref_seq', 'var_seq1', 'var_seq2'])
            common_all = common_all.append(common, ignore_index=True)

        disease_all = disease_all.append(norepeat_tf, ignore_index=True)
        normal_all = normal_all.append(norepeat_nf, ignore_index=True)

    # Function to compute variant presence counts across sample categories
    def compute_subject_counts(df, group_all, common_all, tumorname):
        t0n0, t0n1, t1n0, t1n1 = [], [], [], []
        for i in range(df.shape[0]):
            var = pd.DataFrame(df.iloc[i, :]).T
            t1n1.append(pd.merge(var, common_all, on=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2']).shape[0])
            t1n0.append(pd.merge(var, disease_all, on=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2']).shape[0])
            t0n1.append(pd.merge(var, normal_all, on=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_seq2']).shape[0])
            t0n0.append(len(tumorname) - t1n1[-1] - t1n0[-1] - t0n1[-1])
        df['subt1n1'] = t1n1
        df['subt1n0'] = t1n0
        df['subt0n1'] = t0n1
        df['subt0n0'] = t0n0
        return df

    # Apply subject count function to each variant group
    normal_var = compute_subject_counts(normal_var, normal_all, common_all, tumorname)
    common_var = compute_subject_counts(common_var, normal_all, common_all, tumorname)
    disease_var = compute_subject_counts(disease_var, normal_all, common_all, tumorname)

    # Save results to Excel
    writer = pd.ExcelWriter("./results/all_variants.xlsx")
    disease_var.to_excel(writer, sheet_name="tumor", index=False)
    normal_var.to_excel(writer, sheet_name="normal", index=False)
    common_var.to_excel(writer, sheet_name="common", index=False)
    writer.close()

    print("Excel export complete.")
