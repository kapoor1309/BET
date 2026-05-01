configfile: "config.yaml"

import glob
import os

AA = config["target_aa"]
STRIDE = config.get("stride_path", "/usr/local/bin/stride")

PDBS = [
    os.path.basename(f).replace(".pdb.gz", "")
    for f in glob.glob("pdbs/*.pdb.gz")
]


rule all:
    input:
        "final/ss_profile_HHH_for_arg_with_valid_runs.png",
        "final/valid_pdbs.txt"


rule unzip_pdb:
    input:
        pdb = "pdbs/{pdb}.pdb.gz"
    output:
        unzipped = temp("unzipped_pdbs/{pdb}.pdb")
    shell:
        "zcat {input.pdb} > {output.unzipped}"


rule run_stride:
    input:
        unzipped = "unzipped_pdbs/{pdb}.pdb"
    output:
        ss_out = "stride_out/{pdb}.ss.out"
    params:
        stride = STRIDE
    shell:
        "{params.stride} {input.unzipped} > {output.ss_out} 2>/dev/null || touch {output.ss_out}"


rule extract_context:
    input:
        ss = "stride_out/{pdb}.ss.out"
    output:
        tsv = "contexts/context_for_{aa}_in_{pdb}.tsv"
    params:
        aa = "{aa}"
    shell:
        "python scripts/make_contexts.py {input.ss} {output.tsv} {params.aa}"


rule calculate_angles:
    input:
        expand("contexts/context_for_{aa}_in_{pdb}.tsv", pdb=PDBS, aa=[AA])
    output:
        angles = "final/angles.tsv",
        valid_pdbs = "final/valid_pdbs.txt"
    params:
        aa = AA
    shell:
        "python scripts/compute_angles.py contexts {output.angles} {params.aa} pdbs"


rule plot_angles:
    input:
        "final/angles.tsv"
    output:
        "final/ss_profile_HHH_for_arg_with_valid_runs.png"
    shell:
        "python scripts/make_plot.py {input} {output}"
