import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cygnus as cyg
import glob

sns.set_style('ticks')
sns.set_palette('colorblind')
plt.style.use('../2021.08_plotstyle.mplstyle')
plt.rcParams['font.size'] = 14 # for sbcb seminar

# Plotting as violin plots (apo vs. with protein)
# PLOTTING WT data only here!!!!
# y_data attribute is of interest here
#crystal_dist = [4.62, 4.70, 4.66]

# data in nm so multiply by 10

t24_data = [cyg.XvgFile(file).process_data().y_data*10
            for file in sorted(glob.glob('*t24*WT.xvg'))]
t33a_data = [cyg.XvgFile(file).process_data().y_data*10
             for file in sorted(glob.glob('*t33a*WT.xvg'))]
t33b_data = [cyg.XvgFile(file).process_data().y_data*10
             for file in sorted(glob.glob('*t33b*WT.xvg'))]

protein_labels = ['TRIM24', 'TRIM33\u03b1', 'TRIM33\u03b2']
run_labels = ['Apo', 'Peptide bound']

t24_df = pd.DataFrame(t24_data).T.dropna()
print(t24_df.describe())
t24_df.columns = run_labels
t24_df['Protein'] = protein_labels[0]
t24_df_melted = pd.melt(t24_df, id_vars='Protein')

t33a_df = pd.DataFrame(t33a_data).T.dropna()
print(t33a_df.describe())
t33a_df.columns = run_labels
t33a_df['Protein'] = protein_labels[1]
t33a_df_melted = pd.melt(t33a_df, id_vars='Protein')

t33b_df = pd.DataFrame(t33b_data).T.dropna()
print(t33b_df.describe())
t33b_df.columns = run_labels
t33b_df['Protein'] = protein_labels[2]
t33b_df_melted = pd.melt(t33b_df, id_vars='Protein')

df_list = [t24_df_melted, t33a_df_melted, t33b_df_melted]
combined_df = pd.concat(df_list)
print(combined_df)

#----------------------
# WT violin plot
#----------------------
fig, ax = plt.subplots(figsize=(6,4))
sns.violinplot(x='Protein', y='value', hue='variable',
               data=combined_df, split=True, cut=0)
ax.set(xlabel=None,
       ylabel='Distance (Å)',
       ylim=(3.5,12.5),)
#ax.axhline(y=4.62, color='gray', linestyle='--')
ax.legend(frameon=False, fontsize=12)
sns.despine()
fig.tight_layout()
fig.savefig('WT_violinplot.png', dpi=300)
plt.show()
