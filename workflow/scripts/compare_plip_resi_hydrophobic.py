import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])


def make_df_single_type(file, label, resid=None):
    """For making a DataFrame of PLIP breakdown data (one interaction type only)

    Reads my PLIP full breakdown file and returns a dataframe. This is for the
    hydrophobic contacts by residue sheet, so read sheet_name=3 and skip the first
    column. Uses Protein_ID as index to help with sorting later on.

    Args:
        file (str): Path to excel sheet containing PLIP data.
        label (str): Protein label.
        resid (int): Specify residue number of interest. Defaults to None.

    Returns:
        DataFrame: DataFrame containing per-residue interactions by interaction type.
    """

    df = pd.read_excel(file, sheet_name=3, usecols=[1, 2, 3, 4, 5], index_col=0)
    df.fillna(value=0, inplace=True)
    df.reset_index(inplace=True)
    df.set_index("Protein_ID", inplace=True)
    df["Peptide_number"].astype("int")
    df["Protein"] = label

    if resid is not None:
        return df[df["Peptide_number"] == resid]
    else:
        return df


# comparing T24/T33A/T33B
t24_file = snakemake.input[0]
t33a_file = snakemake.input[1]
t33b_file = snakemake.input[2]
files = [t24_file, t33a_file, t33b_file]

t24_resid = snakemake.params.t24_resid
t33a_resid = snakemake.params.t33a_resid
t33b_resid = snakemake.params.t33b_resid
resids = [t24_resid, t33a_resid, t33b_resid]

peptide_resid = snakemake.params.peptide_resid
labels = snakemake.params.protein_labels

all_df = pd.DataFrame()

for i, (file, resid) in enumerate(zip(files, resids)):
    df = make_df_single_type(file, labels[i], resid=peptide_resid)
    sliced_df = df.loc[df.index.isin(resid)].reindex(resid).reset_index()
    sliced_df["Set"] = pd.factorize(sliced_df["Protein_ID"])[0]
    all_df = pd.concat([all_df, sliced_df])

# changed size for thesis
fig, ax = plt.subplots(figsize=(5, 3.5), constrained_layout=True)
sns.barplot(data=all_df, x="Set", y="size", hue="Protein", ax=ax)

ymin = int(getattr(snakemake.params, "ymin", 0))
ymax = int(getattr(snakemake.params, "ymax", 1000))
# change labels for thesis
xtick_labels = snakemake.params.t33b_labels
"""
xtick_labels = [
    f"{i}/{j}/{k}"
    for (i, j, k) in zip(
        snakemake.params.t24_labels,
        snakemake.params.t33a_labels,
        snakemake.params.t33b_labels,
    )
]
"""
ax.set(xlabel="Protein residue", ylabel="Count", ylim=(ymin, ymax))
ax.set_xticklabels(xtick_labels)
ax.legend(frameon=False, loc="upper left")

fig.savefig(snakemake.output[0], dpi=600)
