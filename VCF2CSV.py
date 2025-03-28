# VCF2CSV.py
# Converts VCF files into separate CSV files for tumor and normal samples.
# Original code by Dayo Shittu.

import os
import pandas as pd
import vcf
# To import vcf, you might need to manually add "PyVCF" in the preferences settings >Project:>Python Interpreter

def writecsv(inputfile):
    # Read VCF file
    vcf_reader = vcf.Reader(open(inputfile, 'r'))

    # Lists to store parsed variant information, one for tumor, one for normal
    normal_var = []
    tumor_var = []

    # Placeholder values used in variant output
    var_score = [28]
    refct = []
    altct = []
    refct1 = []
    altct1 = []

    # To extract patient ID from metadata
    vcfread = vcf_reader.metadata['INDIVIDUAL']
    temp = str(vcfread[0]).split("'")
    inp1 = temp[3].split("-")[2]
    patient_id = [inp1]

    for record in vcf_reader:
        ################################### Tumor sample processing #################################################
        curr_var = []
        left = []
        right = []
        ref_seq = []
        alt_seq = []

        tumorAD = record.samples[1]["AD"]
        normalAD = record.samples[0]["AD"]

        # Check if multiple ALT alleles
        if len(tumorAD) == 3:
            sum_ad = sum(tumorAD)
            x, y, z = tumorAD

            # If at least two alleles have >5% frequency
            if (x / sum_ad) > 0.05 and (y / sum_ad) > 0.05:
                refct.append(x); altct.append(y)
                tumor_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

            elif (x / sum_ad) > 0.05 and (z / sum_ad) > 0.05:
                refct.append(x); altct.append(z)
                tumor_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

            elif (y / sum_ad) > 0.05 and (z / sum_ad) > 0.05:
                refct.append(y); altct.append(z)
                tumor_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

        else:
            # Single ALT allele case
            sum_ad = tumorAD[0] + tumorAD[1]
            refcount, altcount = tumorAD

            if (refcount / sum_ad) > 0.05 and (altcount / sum_ad) > 0.05:
                refct.append(refcount); altct.append(altcount)
                tumor_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

        ########################################## Normal sample processing #########################################
        curr_var1 = []
        left1 = []
        right1 = []
        ref_seq1 = []
        alt_seq1 = []

        if len(normalAD) == 3:
            sum_ad = sum(normalAD)
            x, y, z = normalAD
            if (x / sum_ad) > 0.05 and (y / sum_ad) > 0.05:
                refct1.append(x); altct1.append(y)
                normal_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

            elif (x / sum_ad) > 0.05 and (z / sum_ad) > 0.05:
                refct1.append(x); altct1.append(z)
                normal_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

            elif (y / sum_ad) > 0.05 and (z / sum_ad) > 0.05:
                refct1.append(y); altct1.append(z)
                normal_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

        else:
            sum_ad = normalAD[0] + normalAD[1]
            refcount, altcount = normalAD

            if (refcount / sum_ad) > 0.05 and (altcount / sum_ad) > 0.05:
                refct1.append(refcount); altct1.append(altcount)
                normal_var.append([record.CHROM, record.POS, record.POS + len(record.ALT), record.REF, str(record.ALT[0])] + var_score + patient_id)

    # Convert to DataFrames
    tumor = pd.DataFrame(tumor_var, columns=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_score', 'patient_id'])
    normal = pd.DataFrame(normal_var, columns=['chrom', 'left', 'right', 'ref_seq', 'var_seq1', 'var_score', 'patient_id'])

    # Infer major/minor alleles and add columns
    normal['count1'] = altct1
    normal['count2'] = refct1
    tumor['count1'] = altct
    tumor['count2'] = refct

    # Harmonize var_seq1 and var_seq2 for sorting purposes
    normal['var_seq2'] = [ref if ref > alt else alt for ref, alt in zip(normal['ref_seq'], normal['var_seq1'])]
    tumor['var_seq2'] = [ref if ref > alt else alt for ref, alt in zip(tumor['ref_seq'], tumor['var_seq1'])]

    # Reorder var_seq1 and count1/count2 if swapped
    for df in [normal, tumor]:
        for i in range(df.shape[0]):
            if df.at[i, 'ref_seq'] < df.at[i, 'var_seq1']:
                df.at[i, 'var_seq1'], df.at[i, 'ref_seq'] = df.at[i, 'ref_seq'], df.at[i, 'var_seq1']
                df.at[i, 'count1'], df.at[i, 'count2'] = df.at[i, 'count2'], df.at[i, 'count1']

    # Add variant index
    normal.insert(0, 'var_index', range(normal.shape[0]))
    tumor.insert(0, 'var_index', range(tumor.shape[0]))

    # Save to CSV files
    base = os.path.splitext(os.path.basename(inputfile))[0]
    tumor.to_csv(f"{inp1}_{base}_tumor.csv", index=False)
    normal.to_csv(f"{inp1}_{base}_normal.csv", index=False)

if __name__ == '__main__':
    # Set VCF file directory to read file into defined function 
    os.chdir("./data/vcfs")
    vcf_files = [f for f in os.listdir() if f.endswith('.vcf')]
    for file in vcf_files:
        writecsv(file)
    print("Conversion complete.")
