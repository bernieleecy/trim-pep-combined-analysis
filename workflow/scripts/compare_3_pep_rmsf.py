# Specifically for peptide RMSF

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])


def process_rmsf(file_list):
    """
    For processing a list of xvg files to extract RMSF

    As I want to use df.describe(), don't transpose the data

    Returns the summary statistics of the dataframe (Average is row 1, stdev is
    row 2)
    """
    rmsf = []

    for file in file_list:
        with open(file) as infile:
            val = []

            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    data = line.strip().split("  ")  # two spaces between each entry
                    rmsf_val = float(data[1].strip())
                    val.append(rmsf_val * 10)  # convert to A

            rmsf.append(val)

    df = pd.DataFrame(rmsf).describe()

    return df


files_1 = [snakemake.input.set_1a, snakemake.input.set_1b, snakemake.input.set_1c]
files_2 = [snakemake.input.set_2a, snakemake.input.set_2b, snakemake.input.set_2c]
files_3 = [snakemake.input.set_3a, snakemake.input.set_3b, snakemake.input.set_3c]
resi_types = snakemake.params.types
labels = snakemake.params.protein_labels

ymin = getattr(snakemake.params, "ymin", 0)
ymax = getattr(snakemake.params, "ymax", 5)

# peptide residues hard coded in
resi = np.arange(1, 21)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True, sharey=True)

for i, (set_1, set_2, set_3) in enumerate(zip(files_1, files_2, files_3)):
    # all data for the first set, excludes C terminal cap
    pep_1_rmsf_data = process_rmsf(set_1)
    pep_1_rmsf_avg = pep_1_rmsf_data.iloc[1, :-1]
    pep_1_rmsf_stdev = pep_1_rmsf_data.iloc[2, :-1]

    # all data for the second set, excludes C terminal cap
    pep_2_rmsf_data = process_rmsf(set_2)
    pep_2_rmsf_avg = pep_2_rmsf_data.iloc[1, :-1]
    pep_2_rmsf_stdev = pep_2_rmsf_data.iloc[2, :-1]

    # all data for the third set, excludes C terminal cap
    pep_3_rmsf_data = process_rmsf(set_3)
    pep_3_rmsf_avg = pep_3_rmsf_data.iloc[1, :-1]
    pep_3_rmsf_stdev = pep_3_rmsf_data.iloc[2, :-1]

    # plot on axes
    axes[i].errorbar(
        resi,
        pep_1_rmsf_avg,
        yerr=pep_1_rmsf_stdev,
        alpha=0.7,
        marker="o",
        capsize=3,
        label=resi_types[0],
    )
    axes[i].errorbar(
        resi,
        pep_2_rmsf_avg,
        yerr=pep_2_rmsf_stdev,
        alpha=0.7,
        marker="o",
        capsize=3,
        label=resi_types[1]
    )
    axes[i].errorbar(
        resi,
        pep_3_rmsf_avg,
        yerr=pep_3_rmsf_stdev,
        alpha=0.7,
        marker="o",
        capsize=3,
        label=resi_types[2],
    )

    # updated for thesis
    xticks = [1,5,10,15,20]
    if i == 0:
        axes[i].set(
            xticks=xticks,
            xticklabels=xticks,
            ylim=(ymin, ymax),
        )
        axes[i].set_xlabel("Peptide residue number", fontsize=14)
        axes[i].set_ylabel("RMSF (Å)", fontsize=14)
        axes[i].set_title(labels[i], pad=8, fontsize=16)
        axes[i].legend(frameon=False, loc="upper left", fontsize=13)
    else:
        axes[i].set(
            xticks=xticks,
            xticklabels=xticks,
        )
        axes[i].set_xlabel("Peptide residue number", fontsize=14)
        axes[i].set_title(labels[i], pad=8, fontsize=16)

fig.savefig(snakemake.output[0], dpi=600)
