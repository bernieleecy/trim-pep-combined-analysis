# Snakemake workflow: trim-pep-combined-analysis

A snakemake workflow for analysing combinations of TRIM–peptide simulations.
Most often used for collating T24, T33A and T33B data for the thesis. 

Here, I still use `snake_env`, similar to other workflows.

Unlike the other workflows in the thesis, this repo also contains images!
This is one of the few repos that is small enough to be fully uploaded onto Github, even with figures included. 

## Organisation

A little bit different, as symlinks are now made to the results folders in the T24, T33A and T33B pep analysis workflows.

An example of the organisation is shown below (with the plip folder already present in `combined_plots`). 
```
.
├── README.md
├── combined_plots
│   └── plip
├── config
│   ├── config.yaml
│   ├── samples.csv
│   └── t33a_small_mol_plotstyle.mplstyle
├── data
│   ├── t24 -> ../../t24-pep-analysis/results
│   ├── t33a -> ../../t33a-pep-analysis/results
│   └── t33b -> ../../t33b-pep-analysis/results
└── workflow
    ├── Snakefile
    ├── rules
    └── scripts
```

As of 16/01/2023, `samples.csv` is being used to supply the protein name, but in the future the folder name in the csv file may also be useful. 


# Usage 

As mentioned above, the data folder now contains symlinks to the T24/T33A/T33B results folders.
* In theory, if I update the files in these data folders (e.g. by rerunning the T24/T33A/T33B basic analysis workflows), then Snakemake will also know that it has to rerun the combined analysis

To run the basic analysis:
```
snakemake -np # dry run
snakemake -call
```

Unlike the basic analysis workflows, I will backup the `.png` files generated to github.

# A note on VMD work

For the VMD hydrogen bonds, it was faster to manually select the hydrogen bonds of interest, then use the workflow to plot the heatmaps.
This means that the csv files in `data/vmd_manual_curate` will need to be manually updated if I change the VMD hydrogen bonds code in the earlier analysis workflows.
