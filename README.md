# **SigDyn**

**SigDyn: single-cell mutational signature dynamics**

SigDyn is a workflow designed to study mutational processes and their dynamics in tumors using single-cell RNA-seq data from one timepoint.
It integrates single-cell variant information, phylogenetic reconstruction, and mutational signature analysis to track how mutational processes evolve during tumor evolution.

![SigDyn workflow](SigDyn.png)

---

## **Overview of the workflow**

SigDyn consists of 3 steps:

1. **Single-cell somatic mutation detection**
   Reads are aligned using CellRanger and genotype matrices are inferred using PhylinSic, followed by additional filtering of germline variants and artifacts.

2. **Cancer cell phylogeny inference**
   A phylogenetic tree is reconstructed from variant information using FastTree, followed by (optional) manual rooting to reflect tumor evolution.

3. **Mutational signature dynamics**
   Mutational signatures are assigned to single cells and inferred ancestral nodes, allowing visualization of how signature activities evolve along the tree.

---

## **Requirements**

### **Software**

* `snakemake ≥ 7`
* `conda` (or `mamba`)

### **Automatically managed via Conda**

All Python and tool dependencies are handled automatically using snakemake’s `--use-conda`, including:

* Python ≥ 3.10
* pandas, numpy, ete3
* FastTree
* SigProfilerAssignment
* SigProfilerMatrixGenerator

---

## **Inputs**

### **Required**

### **Reference files**

* `COSMIC_v3.4_SBS_GRCh38.txt` should be downloaded from https://cancer.sanger.ac.uk/signatures/downloads/ and added to the directory resources/

## **Germline variant removal**

SigDyn supports **two strategies** for removing germline variants:

### **Panel of Normals (PON)** *(default in this workflow)*

* Variants present in a population-level panel (e.g. 1000 Genomes) are removed.
* For that, the file should be downloaded from https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON 
* It should then be added to resources/ under the name pon.vcf

### **Matched normal samples (when available)**

* If matched normal single-cell or bulk samples are available, variants homozygous alternate present in normal cells can be filtered instead.
* The filtering script can be adapted by the user to use matched normals rather than a PON.


## **Running the workflow (quick start)**

### **Step 1 — Alignment (using CellRanger) and initial single-cell variant calling**

Starting from FASTQ sample(s), CellRanger should be used for aligment and preprocessing following an appropriate pipeline depending on the input data (see https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis). 
Afterwards, variant calling using the PhylinSic workflow (available at https://github.com/U54Bioinformatics/PhylinSic_Project/tree/main) should be implemented from the first step (_Prepare the reference genome_) until the 8th step (_Call/impute genotypes_). This will allow the user to obtain a "genotypes.txt" that will be the input for the remaining steps of the workflow here described. An example (mock) genotypes.txt file is provided in this repository. 

> A note at the end of this README describes some implementation choices used in our setup of the PhylinSic workflow to improve runtime and ensure compatibility with downstream SigDyn analyses.

### **Step 2 — Varinant filtering, and tree inference**

After having the genotypes.txt file under data/phylinsic/, all under remaining steps that follow use the snakemake pipeline in this repository. 

```bash
snakemake --use-conda --conda-frontend conda --cores 8 stage1_unrooted -p
```

This will:

* Filter germline variants
* Suggest a root candidate cell (highest reference burden)
* Convert genotypes to FASTA
* Infer an unrooted phylogenetic tree

Outputs:

* `results/trees/cell_tree.unrooted`
* `results/trees/root_cell.txt`

---

### **Step 3 — Manual rooting (optional)**

SigDyn **intentionally pauses** here to allow biologically informed rooting.

1. Open the unrooted tree in **FigTree**
2. Root the tree (e.g. using the suggested root cell)
3. Export as (in newick format):

   ```
   results/trees/cell_tree.rooted
   ```

In the case the user does not want to root the tree, they must change the naming of the file "cell_tree.unrooted" to "cell_tree.rooted".

---

### **Step 4 — Resume full analysis**

```bash
snakemake --use-conda --conda-frontend conda --cores 8 stage2_full
```

This will:

* Build mutation profiles for all internal (ancestral) nodes
* Prepare VCF files
* Assign COSMIC mutational signatures (leaves)
* Identify active signatures (present in ≥10% of cells)
* Re-assign signatures using only active signatures (leaves & internal nodes)
* Output final signature exposures for all nodes and single-cells.

---

### **Final note regarding PhylinSic adjustments**
  
For the PhylinSic pipeline, some adjustments were made to both improve reliability of smoothing and efficiency of running the pipeline by changing one rule and paralellizing another. We advise users to do the same.

