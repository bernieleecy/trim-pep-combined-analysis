import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])


def make_df(file, label, start_res=1, end_res=20):
    """For making a DataFrame of PLIP hydrophilic and hydrophobic contacts

    Reads my PLIP hydrophilic and hydrophobic contacts breakdown and processes the data
    for specified residues. Supply the residue numbers as is, do not correct for
    python's zero indexing.

    Args:
        file (str): Path to excel sheet containing PLIP data.
        label (str): Protein label.
        start_res (int): First residue of the peptide. Defaults to 1.
        end_res (int): Last residue of the peptide. Defaults to 20.

    Returns:
        DataFrame: DataFrame containing per-residue hydrophilic-hydrophobic contact
        data.
    """
    df = pd.read_excel(file, sheet_name=0, index_col=0)
    df.fillna(value=0, inplace=True)
    df.reset_index(inplace=True)
    df["Peptide_number"].astype("int")

    # minus one from start_res to handle the zero indexing
    df = df.iloc[start_res - 1 : end_res]
    melt_df = pd.melt(
        df,
        id_vars=["Peptide_number"],
        value_vars=["Hydrophobic", "Hydrophilic"],
        var_name="Contact type",
    )
    melt_df["Protein"] = label
    return melt_df


f1 = snakemake.input.files_1
f2 = snakemake.input.files_2
f3 = snakemake.input.files_3
resi_types = snakemake.params.types
labels = snakemake.params.protein_labels

all_df = pd.DataFrame()

# requires python 3.10 to use strict
for i, files in enumerate([f1, f2, f3]):
    type = resi_types[i]

    for file, label in zip(files, labels, strict=True):
        df = make_df(file, label=label, start_res=1, end_res=20)
        df["Type"] = type
        all_df = pd.concat([all_df, df])

peptide_resid = int(getattr(snakemake.params, "peptide_resid"))

# catplot for interactions
resi_plot = sns.catplot(
    x="Type",
    y="value",
    hue="Contact type",
    col="Protein",
    data=all_df[all_df["Peptide_number"] == peptide_resid],
    order=resi_types,
    hue_order=["Hydrophilic", "Hydrophobic"],
    kind="bar",
    palette=["lightsteelblue", "thistle"],
    height=4.0,
    aspect=0.9,
    legend=False,
)

# legend on leftmost axes only
resi_plot.fig.get_axes()[0].legend(loc="upper left", frameon=False, fontsize=11)

ymin = int(getattr(snakemake.params, "ymin", 0))
ymax = int(getattr(snakemake.params, "ymax", 1000))
resi_plot.fig.get_axes()[0].set(ylim=(ymin, ymax))

(
    resi_plot.set_axis_labels("", "Count")
    .set_titles("{col_name}")
    .set_xticklabels(resi_types)
)

resi_plot.tight_layout()
plt.savefig(snakemake.output[0], dpi=600)
