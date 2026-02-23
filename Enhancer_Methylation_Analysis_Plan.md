# Enhancer Methylation Analysis Plan

**Script**: `identify_hyper_hypo_enhancers_T2T.R`  
**Author**: Cascade AI Assistant  
**Date**: 2026-02-20  

---

## 1. Objective

To analyze enhancer methylation patterns from EM-seq data, integrate with RNA-seq gene expression data, and perform statistical tests to identify significant associations between enhancer methylation and gene regulation.

---

## 2. Data Sources

| Data Type | Source | Description |
|---|---|---|
| **Enhancer Annotations** | `C:/PROJECTS/resource/T2T_CHM13/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T_6col.bed` | BED file with enhancer coordinates |
| **Enhancer-Gene Links** | `C:/PROJECTS/resource/T2T_CHM13/ENCFF912FUA_MCF10A_element_gene_links_thresholded_Engreitz_T2T.bed` | Links enhancers to their target genes |
| **EM-seq DMRs** | `C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/dmrseq` | Differential methylation regions from dmrseq |
| **EM-seq Targeted** | `C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/EMseq/Enhancer_targeted` | Pre-filtered enhancer methylation results |
| **RNA-seq DESeq2** | `C:/PROJECTS/Shane/Harding_250611/T2T_CHM13/RNAseq/deseq2/single_factor_analysis_gene_level` | Differential gene expression results |

---

## 3. Analysis Workflow

### 3.1. Setup and Configuration

1. **Load Libraries**: `data.table`, `GenomicRanges`, `dplyr`, `tidyr`, `rtracklayer`, `openxlsx`
2. **Define Paths**: Set base directories for all input and output data
3. **Set Thresholds**:
   - dmrseq: `pval < 0.05`, `|beta| > 0.5`
   - Targeted: `pvalue < 0.05`, `|meth.diff| > 10`
   - RNA-seq: `padj < 0.05`, `|log2FoldChange| > 1`

### 3.2. Data Loading and Preprocessing

1. **Load Enhancer Annotations**: Read 6-column BED file for coordinates
2. **Load Enhancer-Gene Links**: 
   - Read full BED file with target genes
   - Aggregate by enhancer name to create a mapping of enhancer → comma-separated target genes
3. **Identify Comparisons**: Scan dmrseq directory for comparison folders (containing "_vs_")

### 3.3. Per-Comparison Analysis Loop

For each comparison (e.g., `IR10Gy24h_vs_NIR`):

#### 3.3.1. RNA-seq Integration
1. **Map Comparison Name**: Convert EM-seq names to RNA-seq names (e.g., `IR10Gy24h_vs_NIR` → `10Gy24h_vs_NIR`)
2. **Load DESeq2 Results**: Read `_DESeq2_results_regular.csv`
3. **Identify Regulated Genes**:
   - Filter by `padj < 0.05` and `|log2FoldChange| > 1`
   - Classify as "Up" if `log2FoldChange > 0`, "Down" otherwise
   - Create a lookup map: `gene_name → Regulation`

#### 3.3.2. DMRseq Analysis
1. **Load DMRs**: Read `DMR/DMR_table.tsv`
2. **Filter Significant DMRs**: Apply `pval < 0.05` and `|beta| > 0.5`
3. **Classify Methylation**:
   - **Hypermethylated**: `beta > 0.5`
   - **Hypomethylated**: `beta < -0.5`
4. **Intersect with Enhancers**: Find overlapping enhancer regions

#### 3.3.3. Targeted Analysis
1. **Load Results**: Read `*_enhancer_summary.txt` and `*_promoter_summary.txt`
2. **Filter Significant**: Apply `pvalue < 0.05` and `|meth.diff| > 10`
3. **Classify Methylation**:
   - **Hypermethylated**: `meth.diff > 10`
   - **Hypomethylated**: `meth.diff < -10`
4. **Intersect with Enhancers**: Find overlapping enhancer regions

#### 3.3.4. Union and Annotation
1. **Create Unions**: Combine dmrseq and targeted results for hyper- and hypomethylated enhancers
2. **Merge Annotations**:
   - Add enhancer coordinates
   - Add `TargetGenes` from enhancer-gene links
   - Add `Regulated_TargetGenes` by annotating with RNA-seq data
3. **Count Unique Genes**: Calculate number of unique up/down-regulated genes per enhancer set

#### 3.3.5. Statistical Analysis
1. **Fisher's Exact Test 1**: Association between enhancer methylation and gene regulation
   - **Contingency Table**:
     | | Upregulated | Downregulated |
     |---|---|---|
     | **Hypomethylated** | a | b |
     | **Hypermethylated** | c | d |
2. **Fisher's Exact Test 2**: Association between downregulated genes and hypermethylated enhancers
   - **Contingency Table**:
     | | Downregulated | Not Downregulated |
     |---|---|---|
     | **Hypermethylated** | a | b |
     | **Not Hypermethylated** | c | d |

### 3.4. Output Generation

1. **Per-Comparison Files**:
   - `*_Hypermethylated_Enhancers_Union.csv`
   - `*_Hypomethylated_Enhancers_Union.csv`
   - `*_Fisher_Exact_Methylation_vs_Regulation.txt`
2. **Summary File**:
   - `Enhancer_Methylation_Summary_Counts.csv`: Aggregated counts for all comparisons

---

## 4. Key Functions

### `annotate_rna(target_genes_str, reg_map)`
- **Purpose**: Annotates a comma-separated string of genes with their regulation status
- **Input**: Gene string, regulation map from RNA-seq
- **Output**: Formatted string (e.g., "GeneA(Up), GeneB(Down)")

### `count_reg_genes(enhancer_df)`
- **Purpose**: Counts unique up- and down-regulated genes from an enhancer dataframe
- **Input**: Dataframe with `Regulated_TargetGenes` column
- **Output**: Named vector `c(Up = N, Down = M)`

---

## 5. Output Interpretation

### 5.1. Summary Counts Table
- `Hyper_dmrseq`, `Hypo_dmrseq`: Number of enhancers from dmrseq
- `Hyper_Targeted`, `Hypo_Targeted`: Number of enhancers from targeted analysis
- `Hyper_Union`, `Hypo_Union`: Unique enhancers after union
- `*_Genes_Up`, `*_Genes_Down`: Unique regulated genes per enhancer set

### 5.2. Fisher's Exact Test Results
- **Test 1**: Overall association between methylation direction and gene regulation
- **Test 2**: Specific test for hypermethylation and downregulation
- **Key Metrics**: P-value, Odds Ratio, 95% Confidence Interval

---

## 6. Dependencies and Assumptions

1. **File Structure**: Assumes specific directory structure and naming conventions
2. **Genome Build**: All data must be on the same T2T_CHM13 build
3. **Gene Identifiers**: Consistent gene naming across datasets
4. **Statistical Independence**: Assumes Fisher's exact test is appropriate for the data structure

---

## 7. Potential Extensions

1. **Visualization**: Add volcano plots and heatmaps
2. **Pathway Analysis**: Perform enrichment on regulated gene sets
3. **Comparative Analysis**: Compare across multiple conditions or time points
4. **Machine Learning**: Predict gene regulation from methylation patterns

---

## 8. Quality Control

1. **Verify Overlaps**: Ensure enhancer intersections are correct
2. **Check Counts**: Validate that gene counts are consistent across outputs
3. **Statistical Validity**: Ensure assumptions for Fisher's test are met
4. **Reproducibility**: Document all parameters and random seeds

---

## 9. Execution

```bash
Rscript identify_hyper_hypo_enhancers_T2T.R
```

The script will create all necessary output directories and generate comprehensive results for each comparison.
