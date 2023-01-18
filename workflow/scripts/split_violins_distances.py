import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cygnus as cyg
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

# crystal_dist = [4.62, 4.70, 4.66]


def make_df(data, run_labels, protein_label):
    """Make a DataFrame containing distance data and appropriate labels

    Args:
        data (list of lists): data for the DataFrame
        run_labels (list): list of run types (e.g. ["Apo", "Peptide"])
        protein_label (str): protein used

    Returns:
        DataFrame: Long form DataFrame, identifier var is the protein while run types
        have been unpivoted to the row axis
    """
    df = pd.DataFrame(data).T
    df.columns = run_labels
    df["Protein"] = protein_label
    df_melted = pd.melt(df, id_vars="Protein")

    return df_melted


# 1 set per protein, contains 2 xvg files to compare
# data in nm so multiply by 10
data_1 = [
    cyg.XvgFile(file).process_data().y_data * 10 for file in snakemake.input.set_1
]
data_2 = [
    cyg.XvgFile(file).process_data().y_data * 10 for file in snakemake.input.set_2
]
data_3 = [
    cyg.XvgFile(file).process_data().y_data * 10 for file in snakemake.input.set_3
]

run_types = snakemake.params.types
labels = snakemake.params.protein_labels

all_df = pd.DataFrame()

for i, data in enumerate([data_1, data_2, data_3]):
    df = make_df(data, run_types, labels[i])
    all_df = pd.concat([all_df, df])

ylabel = getattr(snakemake.params, "ylabel", "Distance (Å)")
ymin = getattr(snakemake.params, "ymin", 3)
ymax = getattr(snakemake.params, "ymax", 15)

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
sns.violinplot(x="Protein", y="value", hue="variable", data=all_df, split=True, cut=0)
ax.set(
    xlabel=None,
    ylabel=ylabel,
    ylim=(ymin, ymax),
)
ax.legend(frameon=False, loc="upper left")

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
