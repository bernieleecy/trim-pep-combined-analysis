t24_unmod_rmsf_files = expand(
    "data/t24/unmod_pep/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t24_mod_rmsf_files = expand(
    "data/t24/mod_pep/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t24_mod_trans_rmsf_files = expand(
    "data/t24/mod_pep_trans/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t24_n980a_rmsf_files = expand(
    "data/t24/n980a/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33a_unmod_rmsf_files = expand(
    "data/t33a/unmod_pep/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33a_mod_rmsf_files = expand(
    "data/t33a/mod_pep/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33a_mod_trans_rmsf_files = expand(
    "data/t33a/mod_pep_trans/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33a_n1039a_rmsf_files = expand(
    "data/t33a/n1039a/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33b_unmod_rmsf_files = expand(
    "data/t33b/unmod_pep/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33b_mod_rmsf_files = expand(
    "data/t33b/mod_pep/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33b_mod_trans_rmsf_files = expand(
    "data/t33b/mod_pep_trans/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)
t33b_n1039a_rmsf_files = expand(
    "data/t33b/n1039a/peptide/data/{i}-pep_backbone_rmsf.xvg", i=IDS
)

labels = ["TRIM24", "TRIM33\u03b1", "TRIM33\u03b2"]
ymax_rmsf = 6


rule plot_pep_rmsf_k18_k18ac_cis:
    input:
        set_1a=t24_unmod_rmsf_files,
        set_1b=t33a_unmod_rmsf_files,
        set_1c=t33b_unmod_rmsf_files,
        set_2a=t24_mod_rmsf_files,
        set_2b=t33a_mod_rmsf_files,
        set_2c=t33b_mod_rmsf_files,
    output:
        "combined_plots/peptide/pep_backbone_rmsf_k18_k18ac_cis.png",
    params:
        protein_labels=labels,
        types=["Unmodified peptide", "Modified peptide (C)"],
        ymax=ymax_rmsf,
    script:
        "../scripts/compare_2_pep_rmsf.py"


rule plot_pep_rmsf_k18_k18ac_trans:
    input:
        set_1a=t24_unmod_rmsf_files,
        set_1b=t33a_unmod_rmsf_files,
        set_1c=t33b_unmod_rmsf_files,
        set_2a=t24_mod_trans_rmsf_files,
        set_2b=t33a_mod_trans_rmsf_files,
        set_2c=t33b_mod_trans_rmsf_files,
    output:
        "combined_plots/peptide/pep_backbone_rmsf_k18_k18ac_trans.png",
    params:
        protein_labels=labels,
        types=["Unmodified peptide", "Modified peptide (T)"],
        ymax=ymax_rmsf,
    script:
        "../scripts/compare_2_pep_rmsf.py"


rule plot_pep_rmsf_k18_k18ac_all:
    input:
        set_1a=t24_unmod_rmsf_files,
        set_1b=t33a_unmod_rmsf_files,
        set_1c=t33b_unmod_rmsf_files,
        set_2a=t24_mod_rmsf_files,
        set_2b=t33a_mod_rmsf_files,
        set_2c=t33b_mod_rmsf_files,
        set_3a=t24_mod_trans_rmsf_files,
        set_3b=t33a_mod_trans_rmsf_files,
        set_3c=t33b_mod_trans_rmsf_files,
    output:
        "combined_plots/peptide/pep_backbone_rmsf_k18_k18ac_all.png",
    params:
        protein_labels=labels,
        types=["Unmodified peptide", "Modified peptide (C)", "Modified peptide (T)"],
        ymax=ymax_rmsf,
    script:
        "../scripts/compare_3_pep_rmsf.py"


rule plot_pep_rmsf_k18_k18ac_cis_mut:
    input:
        set_1a=t24_unmod_rmsf_files,
        set_1b=t33a_unmod_rmsf_files,
        set_1c=t33b_unmod_rmsf_files,
        set_2a=t24_mod_rmsf_files,
        set_2b=t33a_mod_rmsf_files,
        set_2c=t33b_mod_rmsf_files,
        set_3a=t24_n980a_rmsf_files,
        set_3b=t33a_n1039a_rmsf_files,
        set_3c=t33b_n1039a_rmsf_files,
    output:
        "combined_plots/peptide/pep_backbone_rmsf_k18_k18ac_cis_mut.png",
    params:
        protein_labels=labels,
        types=[
            "Unmodified peptide",
            "Modified peptide (C)",
            "Modified peptide (C, N to A)",
        ],
        ymax=ymax_rmsf,
    script:
        "../scripts/compare_3_pep_rmsf.py"
