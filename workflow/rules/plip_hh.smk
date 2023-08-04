unmod_files = expand("data/{protein}/unmod_pep/plip/plip_data.xlsx", protein=PROTEINS)
mod_files = expand("data/{protein}/mod_pep/plip/plip_data.xlsx", protein=PROTEINS)
mod_trans_files = expand(
    "data/{protein}/mod_pep_trans/plip/plip_data.xlsx", protein=PROTEINS
)

labels = ["TRIM24", "TRIM33\u03b1", "TRIM33\u03b2"]
ymax_k18 = 1700
ymax_k9 = 1200
ymax_r17 = 900


rule plot_plip_k18_k18ac_cis:
    input:
        files_1=unmod_files,
        files_2=mod_files,
    output:
        "combined_plots/plip_hh/plip_k18_k18ac_cis.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        types=["K18", "K18Ac"],
        ymax=ymax_k18,
    script:
        "../scripts/compare_2_plip_hydrophilic_hydrophobic.py"


rule plot_plip_k18_k18ac_trans:
    input:
        files_1=unmod_files,
        files_2=mod_trans_files,
    output:
        "combined_plots/plip_hh/plip_k18_k18ac_trans.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        types=["K18", "K18Ac"],
        ymax=ymax_k18,
    script:
        "../scripts/compare_2_plip_hydrophilic_hydrophobic.py"


rule plot_plip_k18_k18ac_all:
    input:
        files_1=unmod_files,
        files_2=mod_files,
        files_3=mod_trans_files,
    output:
        "combined_plots/plip_hh/plip_k18_k18ac_all.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        types=["K18", "K18Ac (C)", "K18Ac (T)"],
        ymax=ymax_k18,
    script:
        "../scripts/compare_3_plip_hydrophilic_hydrophobic.py"


rule plot_plip_k9_k9me3_cis:
    """
    Dummy lambda function used to escape the braces around {3}
    """
    input:
        files_1=unmod_files,
        files_2=mod_files,
    output:
        "combined_plots/plip_hh/plip_k9_k9me3_cis.png",
    params:
        peptide_resid=9,
        protein_labels=labels,
        types=lambda wc: ["K9", "K9Me$_\mathregular{3}$"],
        ymax=ymax_k9,
    script:
        "../scripts/compare_2_plip_hydrophilic_hydrophobic.py"


rule plot_plip_k9_k9me3_trans:
    input:
        files_1=unmod_files,
        files_2=mod_trans_files,
    output:
        "combined_plots/plip_hh/plip_k9_k9me3_trans.png",
    params:
        peptide_resid=9,
        protein_labels=labels,
        types=lambda wc: ["K9", "K9Me$_\mathregular{3}$"],
        ymax=ymax_k9,
    script:
        "../scripts/compare_2_plip_hydrophilic_hydrophobic.py"


rule plot_plip_k9_k9me3_all:
    input:
        files_1=unmod_files,
        files_2=mod_files,
        files_3=mod_trans_files,
    output:
        "combined_plots/plip_hh/plip_k9_k9me3_all.png",
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
        "../scripts/compare_3_plip_hydrophilic_hydrophobic.py"


rule plot_plip_r17_cis:
    input:
        files_1=unmod_files,
        files_2=mod_files,
    output:
        "combined_plots/plip_hh/plip_r17_cis.png",
    params:
        peptide_resid=17,
        protein_labels=labels,
        types=["K18", "K18Ac"],
        ymax=ymax_r17,
    script:
        "../scripts/compare_2_plip_hydrophilic_hydrophobic.py"


rule plot_plip_r17_trans:
    input:
        files_1=unmod_files,
        files_2=mod_trans_files,
    output:
        "combined_plots/plip_hh/plip_r17_trans.png",
    params:
        peptide_resid=17,
        protein_labels=labels,
        types=["K18", "K18Ac"],
        ymax=ymax_r17,
    script:
        "../scripts/compare_2_plip_hydrophilic_hydrophobic.py"


rule plot_plip_r17_all:
    input:
        files_1=unmod_files,
        files_2=mod_files,
        files_3=mod_trans_files,
    output:
        "combined_plots/plip_hh/plip_r17_all.png",
    params:
        peptide_resid=17,
        protein_labels=labels,
        types=["K18", "K18Ac (C)", "K18Ac (T)"],
        ymax=ymax_r17,
    script:
        "../scripts/compare_3_plip_hydrophilic_hydrophobic.py"
