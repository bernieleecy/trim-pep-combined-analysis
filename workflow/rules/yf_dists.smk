# T24 files
t24_apo_f979_file = "data/t24/apo/f979_dists/data/f979_y935.xvg"
t24_unmod_f979_file = "data/t24/unmod_pep/f979_dists/data/f979_y935.xvg"
t24_mod_f979_file = "data/t24/mod_pep/f979_dists/data/f979_y935.xvg"
t24_mod_trans_f979_file = "data/t24/mod_pep_trans/f979_dists/data/f979_y935.xvg"
t24_n980a_f979_file = "data/t24/n980a/f979_dists/data/f979_y935.xvg"
# T33A files
t33a_apo_f1038_file = "data/t33a/apo/f1038_dists/data/f1038_y993.xvg"
t33a_unmod_f1038_file = "data/t33a/unmod_pep/f1038_dists/data/f1038_y993.xvg"
t33a_mod_f1038_file = "data/t33a/mod_pep/f1038_dists/data/f1038_y993.xvg"
t33a_mod_trans_f1038_file = "data/t33a/mod_pep_trans/f1038_dists/data/f1038_y993.xvg"
t33a_n1039a_f1038_file = "data/t33a/n1039a/f1038_dists/data/f1038_y993.xvg"
# T33B files
t33b_apo_f1038_file = "data/t33b/apo/f1038_dists/data/f1038_y993.xvg"
t33b_unmod_f1038_file = "data/t33b/unmod_pep/f1038_dists/data/f1038_y993.xvg"
t33b_mod_f1038_file = "data/t33b/mod_pep/f1038_dists/data/f1038_y993.xvg"
t33b_mod_trans_f1038_file = "data/t33b/mod_pep_trans/f1038_dists/data/f1038_y993.xvg"
t33b_n1039a_f1038_file = "data/t33b/n1039a/f1038_dists/data/f1038_y993.xvg"

labels = ["TRIM24", "TRIM33\u03b1", "TRIM33\u03b2"]
ymin_yf = 3.5
ymax_yf = 12


rule plot_yf_apo_k18ac_cis:
    input:
        set_1=[t24_apo_f979_file, t24_mod_f979_file],
        set_2=[t33a_apo_f1038_file, t33a_mod_f1038_file],
        set_3=[t33b_apo_f1038_file, t33b_mod_f1038_file],
    output:
        "combined_plots/YF_dists/YF_apo_k18ac_cis_violins.png",
    params:
        protein_labels=labels,
        types=["Apo", "Peptide bound (C)"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"


rule plot_yf_apo_k18ac_trans:
    input:
        set_1=[t24_apo_f979_file, t24_mod_trans_f979_file],
        set_2=[t33a_apo_f1038_file, t33a_mod_trans_f1038_file],
        set_3=[t33b_apo_f1038_file, t33b_mod_trans_f1038_file],
    output:
        "combined_plots/YF_dists/YF_apo_k18ac_trans_violins.png",
    params:
        protein_labels=labels,
        types=["Apo", "Peptide bound (T)"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"


rule plot_yf_apo_k18ac_mut:
    input:
        set_1=[t24_apo_f979_file, t24_n980a_f979_file],
        set_2=[t33a_apo_f1038_file, t33a_n1039a_f1038_file],
        set_3=[t33b_apo_f1038_file, t33b_n1039a_f1038_file],
    output:
        "combined_plots/YF_dists/YF_apo_k18ac_mut_violins.png",
    params:
        protein_labels=labels,
        types=["Apo", "Peptide bound (C), N to A"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"


rule plot_yf_unmod_k18ac_cis:
    input:
        set_1=[t24_unmod_f979_file, t24_mod_f979_file],
        set_2=[t33a_unmod_f1038_file, t33a_mod_f1038_file],
        set_3=[t33b_unmod_f1038_file, t33b_mod_f1038_file],
    output:
        "combined_plots/YF_dists/YF_unmod_k18ac_cis_violins.png",
    params:
        protein_labels=labels,
        types=["Peptide bound (unmod)", "Peptide bound (C)"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"


rule plot_yf_unmod_k18ac_trans:
    input:
        set_1=[t24_unmod_f979_file, t24_mod_trans_f979_file],
        set_2=[t33a_unmod_f1038_file, t33a_mod_trans_f1038_file],
        set_3=[t33b_unmod_f1038_file, t33b_mod_trans_f1038_file],
    output:
        "combined_plots/YF_dists/YF_unmod_k18ac_trans_violins.png",
    params:
        protein_labels=labels,
        types=["Peptide bound (unmod)", "Peptide bound (T)"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"


rule plot_yf_unmod_k18ac_mut:
    input:
        set_1=[t24_unmod_f979_file, t24_n980a_f979_file],
        set_2=[t33a_unmod_f1038_file, t33a_n1039a_f1038_file],
        set_3=[t33b_unmod_f1038_file, t33b_n1039a_f1038_file],
    output:
        "combined_plots/YF_dists/YF_unmod_k18ac_mut_violins.png",
    params:
        protein_labels=labels,
        types=["Peptide bound (unmod)", "Peptide bound (C), N to A"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"


rule plot_yf_k18ac_cis_trans:
    input:
        set_1=[t24_mod_f979_file, t24_mod_trans_f979_file],
        set_2=[t33a_mod_f1038_file, t33a_mod_trans_f1038_file],
        set_3=[t33b_mod_f1038_file, t33b_mod_trans_f1038_file],
    output:
        "combined_plots/YF_dists/YF_k18ac_cis_trans_violins.png",
    params:
        protein_labels=labels,
        types=["Peptide bound (C)", "Peptide bound (T)"],
        ymin=ymin_yf,
        ymax=ymax_yf,
    script:
        "../scripts/split_violins_distances.py"
