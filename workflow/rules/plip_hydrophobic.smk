unmod_files = [f"data/{protein}/unmod_pep/plip/plip_data.xlsx" for protein in PROTEINS]
mod_files = [f"data/{protein}/mod_pep/plip/plip_data.xlsx" for protein in PROTEINS]
mod_trans_files = [
    f"data/{protein}/mod_pep_trans/plip/plip_data.xlsx" for protein in PROTEINS
]

labels = ["TRIM24", "TRIM33\u03b1", "TRIM33\u03b2"]
t24_hydrophobic_3 = ["PHE924", "VAL928", "VAL932", "PHE979", "VAL986"]
t24_hydrophobic_1 = ["F924", "V928", "V932", "F979", "V986"]
t33a_hydrophobic_3 = ["PHE982", "VAL986", "ILE990", "PHE1038", "VAL1062"]
t33a_hydrophobic_1 = ["F982", "V986", "I990", "F1038", "V1062"]
t33b_hydrophobic_3 = ["PHE982", "VAL986", "ILE990", "PHE1038", "VAL1045"]
t33b_hydrophobic_1 = ["F982", "V986", "I990", "F1038", "V1045"]

ymax_k18 = 600


rule plot_plip_k18_cis_key_hydrophobic:
    input:
        mod_files,
    output:
        "combined_plots/plip_hydrophobic/plip_k18ac_cis_key_hydrophobic.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        t24_resid=t24_hydrophobic_3,
        t33a_resid=t33a_hydrophobic_3,
        t33b_resid=t33b_hydrophobic_3,
        t24_labels=t24_hydrophobic_1,
        t33a_labels=t33a_hydrophobic_1,
        t33b_labels=t33b_hydrophobic_1,
        ymax=ymax_k18,
    script:
        "../scripts/compare_plip_resi_hydrophobic.py"


rule plot_plip_k18_trans_key_hydrophobic:
    input:
        mod_trans_files,
    output:
        "combined_plots/plip_hydrophobic/plip_k18ac_trans_key_hydrophobic.png",
    params:
        peptide_resid=18,
        protein_labels=labels,
        t24_resid=t24_hydrophobic_3,
        t33a_resid=t33a_hydrophobic_3,
        t33b_resid=t33b_hydrophobic_3,
        t24_labels=t24_hydrophobic_1,
        t33a_labels=t33a_hydrophobic_1,
        t33b_labels=t33b_hydrophobic_1,
        ymax=ymax_k18,
    script:
        "../scripts/compare_plip_resi_hydrophobic.py"
