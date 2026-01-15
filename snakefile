########################################
# Config
########################################

configfile: "config/config.yaml"

PHYLINSIC_GENOTYPES = config["inputs"]["phylinsic_genotypes"]

########################################
# Final targets
########################################

rule all:
    input:
        "results/trees/cell_tree.rooted.named",
        "results/signatures/exposures_active_signatures.txt"

########################################
# Step 1 (conceptual): PhylinSic (external)
########################################

rule phylinsic_external:
    output:
        PHYLINSIC_GENOTYPES
    message:
        """
        ============================================================
        PhylinSic is an external prerequisite for SigDyn.

        Please run PhylinSic separately and place its output at:

            {output}

        See the README for details on required format and references.
        ============================================================
        """
    shell:
        "exit 1"

########################################
# Step 1b: Additional variant filtering
########################################

rule filter_variants:
    input:
        PHYLINSIC_GENOTYPES
    output:
        "results/genotypes/genotypes_no_PON.txt"
    log:
        "results/logs/filter_variants.log"
    conda:
        "envs/base.yaml"
    script:
        "scripts/filter_variants.py"

rule suggest_root_cell:
    input:
        "results/genotypes/genotypes_no_PON.txt"
    output:
        "results/trees/root_cell.txt"
    conda:
        "envs/base.yaml"
    script:
        "scripts/suggest_root_cell.py"


########################################
# Step 2: Phylogeny inference
########################################

rule genotype_to_fasta:
    input:
        "results/genotypes/genotypes_no_PON.txt"
    output:
        "results/trees/mutations.fa"
    conda:
        "envs/base.yaml"
    script:
        "scripts/genotype_to_fasta.py"

rule infer_tree:
    input:
        "results/trees/mutations.fa"
    output:
        "results/trees/cell_tree.unrooted"
    conda:
        "envs/fasttree.yaml"
    shell:
        """
        FastTree -nt {input} > {output}
        """


########################################
# Manual rooting checkpoint
########################################

rule require_rooted_tree:
    input:
        unrooted="results/trees/cell_tree.unrooted",
        root_tip="results/trees/root_cell.txt"
    output:
        touch("results/trees/.rooting_done")
    priority: 100
    run:
        import os

        rooted = "results/trees/cell_tree.rooted"

        # read suggested tip (if file exists)
        suggested = "UNKNOWN"
        if os.path.exists(input.root_tip):
            with open(input.root_tip) as f:
                suggested = f.read().strip()

        if not os.path.exists(rooted):
            msg = f"""
============================================================
MANUAL ROOTING REQUIRED

Unrooted tree:
  {input.unrooted}

Suggested root tip (max RR):
  {suggested}
(full suggestion file: {input.root_tip})

ACTION (FigTree):
  1) Open: {input.unrooted}
  2) Select the tip named: {suggested}
  3) Root the tree on that tip (or per manuscript rule)
  4) Export the tree as: {rooted} in newick format

Then rerun:
    snakemake --use-conda --conda-frontend conda --cores 8 stage2_full
============================================================
"""
            print(msg, flush=True)
            raise ValueError(msg)

        # rooted exists → create sentinel and continue
        open(output[0], "w").close()




########################################
# Step 3: Ancestral mutation profiles
########################################

rule build_ancestral_profiles:
    input:
        rooting_done="results/trees/.rooting_done",
        genotypes="results/genotypes/genotypes_no_PON.txt"
    output:
        genotypes_all_nodes="results/genotypes/genotypes_all_nodes.txt",
        tree_named="results/trees/cell_tree.rooted.named"
    log:
        "results/logs/build_ancestral_profiles.log"
    conda:
        "envs/base.yaml"
    script:
        "scripts/build_ancestral_profiles.py"


########################################
# Step 4: Mutational signature dynamics
########################################

rule prepare_vcfs_leaves:
    input:
        genotypes="results/genotypes/genotypes_no_PON.txt"
    output:
        global_dir=directory("results/vcf_sig_leaves"),
        per_cell_dir=directory("results/vcf_per_cell_leaves")
    conda:
        "envs/base.yaml"
    script:
        "scripts/prepare_vcfs_leaves.py"

rule install_sigprofiler_genome_grch38:
    output:
        touch("resources/sigprofiler_genomes/GRCh38.installed")
    conda:
        "envs/sigprofiler.yaml"
    script:
        "scripts/install_sigprofiler_genome.py"


rule assign_signatures_leaves:
    input:
        genome_installed="resources/sigprofiler_genomes/GRCh38.installed",
        global_vcf_dir="results/vcf_sig_leaves",
        per_cell_vcf_dir="results/vcf_per_cell_leaves"
    output:
        global_out=directory("results/signatures/all_leaves"),
        per_cell_out=directory("results/signatures/per_leaf")
    log:
        "results/logs/assign_signatures_leaves.log"
    params:
        genome_build="GRCh38",
        context="96",
        cosmic_version="3.4"
    conda:
        "envs/sigprofiler.yaml"
    script:
        "scripts/assign_signatures_leaves.py"


rule prepare_vcfs_all_nodes:
    input:
        genotypes="results/genotypes/genotypes_all_nodes.txt"
    output:
        global_dir=directory("results/vcf_sig_all"),
        per_cell_dir=directory("results/vcf_per_cell_all")
    conda:
        "envs/base.yaml"
    script:
        "scripts/prepare_vcfs_leaves.py"


rule make_custom_cosmic_signatures:
    input:
        exposures_dir="results/signatures/per_leaf",
        cosmic="resources/COSMIC_v3.4_SBS_GRCh38.txt"
    output:
        "results/signatures/Custom_COSMIC_SBS_v3.4.txt"
    log:
        "results/logs/make_custom_cosmic.log"
    params:
        min_fraction=0.10
    conda:
        "envs/base.yaml"
    script:
        "scripts/make_custom_cosmic.py"




rule assign_active_signatures:
    input:
        vcf_dir="results/vcf_per_cell_all",
        custom_cosmic="results/signatures/Custom_COSMIC_SBS_v3.4.txt"
    output:
        "results/signatures/exposures_active_signatures.txt"
    params:
        genome="GRCh38"
    conda:
        "envs/sigprofiler.yaml"
    script:
        "scripts/assign_active_signatures.py"

########################################
# Convenience targets
########################################

rule stage1_unrooted:
    input:
        "results/trees/cell_tree.unrooted",
        "results/trees/root_cell.txt"


rule stage2_full:
    input:
        "results/trees/.rooting_done",
        "results/trees/cell_tree.rooted.named",
        "results/signatures/exposures_active_signatures.txt"

rule stage2_rooting_instructions:
    input:
        "results/trees/.rooting_done"
