import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

labels = snakemake.params.protein_labels
df = pd.read_csv(snakemake.input[0], index_col=0)
df.columns = labels

fig_height = df.shape[0] * 1.25

if fig_height <= 10:
    fig, ax = plt.subplots(figsize=(5,fig_height), constrained_layout=True)
else:
    fig, ax = plt.subplots(figsize=(5,10), constrained_layout=True)

ax = sns.heatmap(df, vmin=0, vmax=100,
                 cmap='crest', annot=True, linewidths=0.2,
                 fmt=".1f",
                 cbar_kws={'label': 'Occupancy (%)'},
                 annot_kws={"size": 16})
ax.set(ylabel='')
ax.tick_params(axis='y', rotation=0)
ax.xaxis.set_ticks_position('top')
ax.tick_params(top=False, left=False)
ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=16, fontweight='bold')
ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=16)

fig.savefig(snakemake.output[0], dpi=600)
