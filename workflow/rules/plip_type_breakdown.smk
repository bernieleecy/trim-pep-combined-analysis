unmod_files = expand("data/{protein}/unmod_pep/plip/plip_data.xlsx", protein=PROTEINS)
mod_files = expand("data/{protein}/mod_pep/plip/plip_data.xlsx", protein=PROTEINS)
mod_trans_files = expand(
    "data/{protein}/mod_pep_trans/plip/plip_data.xlsx", protein=PROTEINS
)

labels = ["TRIM24", "TRIM33\u03b1", "TRIM33\u03b2"]
ymax_k18 = 1700
ymax_k9 = 700
ymax_r17 = 600


rule plot_plip_breakdown_k18_k18ac_cis:
    input:
        files_1=unmod_files,
        files_2=mod_files,
    output:
        "combined_plots/plip_breakdown/plip_k18_k18ac_cis_breakdown.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        types=["K18", "K18Ac"],
        ymax=ymax_k18,
    script:
        "../scripts/compare_2_plip_breakdown.py"


rule plot_plip_breakdown_k18_k18ac_trans:
    input:
        files_1=unmod_files,
        files_2=mod_trans_files,
    output:
        "combined_plots/plip_breakdown/plip_k18_k18ac_trans_breakdown.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        types=["K18", "K18Ac"],
        ymax=ymax_k18,
    script:
        "../scripts/compare_2_plip_breakdown.py"


rule plot_plip_breakdown_k18_k18ac_all:
    input:
        files_1=unmod_files,
        files_2=mod_files,
        files_3=mod_trans_files,
    output:
        "combined_plots/plip_breakdown/plip_k18_k18ac_all_breakdown.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        types=["K18", "K18Ac (C)", "K18Ac (T)"],
        ymax=ymax_k18,
    script:
        "../scripts/compare_3_plip_breakdown.py"


rule plot_plip_breakdown_k9_k9me3_cis:
    """
    Dummy lambda function used to escape the braces around {3}
    """
    input:
        files_1=unmod_files,
        files_2=mod_files,
    output:
        "combined_plots/plip_breakdown/plip_k9_k9me3_cis_breakdown.png",
    params:
        peptide_resid=9,
        protein_labels=labels,
        types=lambda wc: ["K9", "K9Me$_\mathregular{3}$"],
        ymax=ymax_k9,
    script:
        "../scripts/compare_2_plip_breakdown.py"


rule plot_plip_breakdown_k9_k9me3_trans:
    input:
        files_1=unmod_files,
        files_2=mod_trans_files,
    output:
        "combined_plots/plip_breakdown/plip_k9_k9me3_trans_breakdown.png",
    params:
        peptide_resid=9,
        protein_labels=labels,
        types=lambda wc: ["K9", "K9Me$_\mathregular{3}$"],
        ymax=ymax_k9,
    script:
        "../scripts/compare_2_plip_breakdown.py"


rule plot_plip_breakdown_k9_k9me3_all:
    input:
        files_1=unmod_files,
        files_2=mod_files,
        files_3=mod_trans_files,
    output:
        "combined_plots/plip_breakdown/plip_k9_k9me3_all_breakdown.png",
    params:
        peptide_resid=9,
        protein_labels=labels,
        types=lambda wc: [
            "K9",
            "K9Me$_\mathregular{3}$ (C)",
            "K9Me$_\mathregular{3}$ (T)",
        ],
        ymax=ymax_k9,
    script:
        "../scripts/compare_3_plip_breakdown.py"


rule plot_plip_breakdown_r17_cis:
    input:
        files_1=unmod_files,
        files_2=mod_files,
    output:
        "combined_plots/plip_breakdown/plip_r17_cis_breakdown.png",
    params:
        peptide_resid=17,
        protein_labels=labels,
        types=["Unmod", "Mod"],
        ymax=ymax_r17,
    script:
        "../scripts/compare_2_plip_breakdown.py"


rule plot_plip_breakdown_r17_trans:
    input:
        files_1=unmod_files,
        files_2=mod_trans_files,
    output:
        "combined_plots/plip_breakdown/plip_r17_trans_breakdown.png",
    params:
        peptide_resid=17,
        protein_labels=labels,
        types=["Unmod", "Mod"],
        ymax=ymax_r17,
    script:
        "../scripts/compare_2_plip_breakdown.py"


rule plot_plip_breakdown_r17_all:
    input:
        files_1=unmod_files,
        files_2=mod_files,
        files_3=mod_trans_files,
    output:
        "combined_plots/plip_breakdown/plip_r17_all_breakdown.png",
    params:
        peptide_resid=17,
        protein_labels=labels,
        types=["Unmod", "Mod (C)", "Mod (T)"],
        ymax=ymax_r17,
    script:
        "../scripts/compare_3_plip_breakdown.py"
