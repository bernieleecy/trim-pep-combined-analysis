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
resi_types = snakemake.params.types
labels = snakemake.params.protein_labels

ymin = getattr(snakemake.params, "ymin", 0)
ymax = getattr(snakemake.params, "ymax", 5)

# peptide residues hard coded in
resi = np.arange(1, 21)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True, sharey=True)

for i, (set_1, set_2) in enumerate(zip(files_1, files_2)):
    # all data for the first set, excludes C terminal cap
    rmsf_data = process_rmsf(set_1)
    pep_1_rmsf_avg = rmsf_data.iloc[1, :-1]
    pep_1_rmsf_stdev = rmsf_data.iloc[2, :-1]

    # all data for the second set, excludes C terminal cap
    mod_pep_data = process_rmsf(set_2)
    pep_2_rmsf_avg = mod_pep_data.iloc[1, :-1]
    pep_2_rmsf_stdev = mod_pep_data.iloc[2, :-1]

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
        label=resi_types[1],
    )

    if i == 0:
        axes[i].set(
            xlabel="Peptide residue number",
            ylabel="RMSF (Å)",
            xticks=resi,
            xticklabels=resi,
            ylim=(ymin, ymax),
            title=labels[i],
        )
        axes[i].tick_params(axis="x", which="major", labelsize=10)
        axes[i].legend(frameon=False, loc="upper left")
    else:
        axes[i].set(
            xlabel="Peptide residue number",
            xticks=resi,
            xticklabels=resi,
            title=labels[i],
        )
        axes[i].tick_params(axis="x", which="major", labelsize=10)

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
