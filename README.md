# **SigDyn**

**SigDyn: single-cell mutational signature dynamics**

SigDyn is a workflow designed to study mutational processes and their dynamics in tumors using single-cell RNA-seq data from one timepoint.
It integrates single-cell variant information, phylogenetic reconstruction, and mutational signature analysis to track how mutational processes evolve during tumor evolution.

![SigDyn workflow](SigDyn.png)

---

## **Overview of the workflow**

SigDyn consists of **three conceptual steps**:

1. **Single-cell somatic mutation detection**
   Genotype matrices are inferred externally using **PhylinSic**, followed by additional filtering of germline variants and artifacts.

2. **Cancer cell phylogeny inference**
   A phylogenetic tree is reconstructed from variant information using FastTree, followed by manual rooting to reflect tumor evolution.

3. **Mutational signature dynamics**
   Mutational signatures are assigned to single cells and inferred ancestral nodes, allowing visualization of how signature activities evolve along the tree.

---

## **Requirements**

### **Software**

* `snakemake ≥ 7`
* `conda` (or `mamba`)

### **Automatically managed via Conda**

All Python and tool dependencies are handled automatically using Snakemake’s `--use-conda`, including:

* Python ≥ 3.10
* pandas, numpy, ete3
* FastTree
* SigProfilerAssignment
* SigProfilerMatrixGenerator

---

## **Inputs**

### **Required**

* **PhylinSic genotype matrix** (`genotypes.txt`)
  A tab-separated file with:

  * rows = variant sites (`chr_pos_ref_alt`)
  * columns = cells
  * entries = `R` (for homozygous reference sites), `H` (for heterozygous sites), or `A` (for homozygous alternate sites)

> The user must run PhylinSic (also available as a Snakemake workflow here: https://github.com/U54Bioinformatics/PhylinSic_Project/tree/main, and provide the genotypes.txt output under data/phylinsic/.
> A note at the end of this README describes some implementation choices used in our setup of the PhylinSic workflow to improve runtime and ensure compatibility with downstream SigDyn analyses.

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


## **Running the workflow**

### **Step 1 — Generate the genotypes.txt file using the snakemake workflow of PhylinSic**


### **Step 2 — Generate unrooted tree & rooting suggestion**

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

### **Step 3 — Manual rooting (required)**

SigDyn **intentionally pauses** here to allow biologically informed rooting.

1. Open the unrooted tree in **FigTree**
2. Root the tree (e.g. using the suggested root cell)
3. Export as (in newick format):

   ```
   results/trees/cell_tree.rooted
   ```

This step ensures that the root reflects tumor evolutionary history rather than an arbitrary midpoint.

---

### **Step 4 — Resume full analysis**

```bash
snakemake --use-conda --conda-frontend conda --cores 8 stage2_full
```

This will:

* Build mutation profiles for **all internal (ancestral) nodes**
* Prepare VCF files
* Assign COSMIC mutational signatures
* Identify **active signatures** (present in ≥10% of cells)
* Re-assign signatures using only active signatures
* Output final signature exposures for all nodes and single-cells.


---

### **Final note regarding PhylinSic adjustments**

