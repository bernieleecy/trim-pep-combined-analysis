import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])


def make_df_breakdown(file, label):
    """For making a DataFrame of PLIP breakdown data

    Reads my PLIP full breakdown file and returns a dataframe grouped by peptide_number
    and type. Reads sheet 2 (sheet_name=1) and skips column 0.

    Args:
        file (str): Path to excel sheet containing PLIP data.
        label (str): Protein label.

    Returns:
        DataFrame: DataFrame containing per-residue interactions by interaction type.
    """

    df = pd.read_excel(file, sheet_name=1, usecols=[1, 2, 3, 4, 5], index_col=0)
    df.fillna(value=0, inplace=True)
    df.reset_index(inplace=True)
    df["Peptide_number"].astype("int")

    grouped_df = df.groupby(by=["Peptide_number", "Type"]).sum()
    grouped_df["Protein"] = label

    return grouped_df


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
        df = make_df_breakdown(file, label=label)
        df["Run type"] = type
        all_df = pd.concat([all_df, df])

# restore Peptide_number and Type as columns
all_df.reset_index(inplace=True)

hue_order = [
    "hydrogen_bonds",
    "salt_bridges",
    "pi_cation_interactions",
    "hydrophobic_interactions",
]
peptide_resid = int(getattr(snakemake.params, "peptide_resid"))

# catplot for interactions
resi_plot = sns.catplot(
    x="Run type",
    y="size",
    hue="Type",
    col="Protein",
    data=all_df[all_df["Peptide_number"] == peptide_resid],
    order=resi_types,
    hue_order=hue_order,
    kind="bar",
    ci=None,
    height=4.0,
    aspect=0.9,
    legend=False,
    facet_kws=dict(despine=False),
)

# legend on leftmost axes only
leg_labels = ["Hydrogen bond", "Salt bridge", "Cation\u2013pi", "Hydrophobic"]
resi_plot.fig.get_axes()[0].legend(
    loc="upper left", frameon=False, fontsize=11, labels=leg_labels
)

ymin = int(getattr(snakemake.params, "ymin", 0))
ymax = int(getattr(snakemake.params, "ymax", 1000))
resi_plot.fig.get_axes()[0].set(ylim=(ymin, ymax))

(
    resi_plot.set_axis_labels("", "Count")
    .set_titles("{col_name}", pad=8)
    .set_xticklabels(resi_types)
)

resi_plot.tight_layout()
plt.savefig(snakemake.output[0], dpi=600)
