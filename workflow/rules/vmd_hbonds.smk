labels = ["TRIM24", "TRIM33\u03b1", "TRIM33\u03b2"]

rule plot_k18ac_vmd_heatmaps:
    input:
        "data/vmd_manual_curate/k18ac_{ct}_{resi}_hbonds.csv"
    output:
        "combined_plots/vmd_hbonds/k18ac_{ct}_{resi}_heatmap.png"
    params:
        protein_labels=labels,
    script:
        "../scripts/vmd_hbonds_heatmap.py"
