#README

Gene Expression Analysis from 

_Environmental memory gained from exposure to extreme pCO2 variability promotes coral cellular acid–base homeostasis_
Kristen T. Brown, Matheus A. Mello-Athayde, Eugenia M. Sampayo, Aaron Chai, Sophie Dove and Katie L. Barott
Published:14 September 2022 [manuscript link](https://doi.org/10.1098/rspb.2022.0941)

_"Ocean acidification is a growing threat to coral growth and the accretion of coral reef ecosystems. Corals inhabiting environments that already endure extreme diel pCO2 fluctuations, however, may represent acidification-resilient populations capable of persisting on future reefs. Here, we examined the impact of pCO2 variability on the reef-building coral Pocillopora damicornis originating from reefs with contrasting environmental histories (variable reef flat versus stable reef slope) following reciprocal exposure to stable (218 ± 9) or variable (911 ± 31) diel pCO2 amplitude (μtam) in aquaria over eight weeks. Endosymbiont density, photosynthesis and net calcification rates differed between origins but not treatment, whereas primary calcification (extension) was affected by both origin and acclimatization to novel pCO2 conditions. At the cellular level, corals from the variable reef flat exhibited less intracellular pH (pHi) acidosis and faster pHi recovery rates in response to experimental acidification stress (pH 7.40) than corals originating from the stable reef slope, suggesting environmental memory gained from lifelong exposure to pCO2 variability led to an improved ability to regulate acid–base homeostasis. These results highlight the role of cellular processes in maintaining acidification resilience and suggest that prior exposure to pCO2 variability may promote more acidification-resilient coral populations in a changing climate."_

## Data

**Sample metadata**
- [Metadata (Origin/Treatment) for all samples](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/TagSeq_Submission/RNA%20Submission%20Sample%20List%20metadata.csv)
- [Physiological Data](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/data/Heron%20pHi%20coral%20physiology%20and%20respirometry%20R_reduced.csv) from [Brown et al., 2022](https://royalsocietypublishing.org/doi/10.1098/rspb.2022.0941)

**Summary of all laboratory work for DNA/RNA extraction:**
- [DNA/RNA Extractions](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/Project-Summary-Barott-and-Brown-Pdam-RNA-DNA-Extractions.md)

**TagSeq Data QC and Processing**
- [Bioinformatic pipeline](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/Heron-Pdam-gene-expression.md)   
- QC of [Raw](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/data/raw_qc) and [Trimmed](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/data/trimmed_qc) reads
- [Raw data: Sequence upload to SRA, BioProject PRJNA934298](https://www.ncbi.nlm.nih.gov/sra/PRJNA934298)
- [Genome mapping rate per sample](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/TagSeq_Output/mapped_reads_counts_Pacuta.txt)
- [Gene count matrix](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/TagSeq_Output/HeronPdam_gene_count_matrix.csv)

**WGCNA analysis**
- [WGCNA Script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/WGCNA/WGCNA.Rmd)
- [GO Enrichment Script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/WGCNA/GO%20analysis.Rmd)
- [Outputs](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/WGCNA)
  - [GO Enrichment](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/WGCNA/GO_analysis)
  - [**Figure 3**](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/WGCNA/Both_with%20phys%20and%20pHi_heatmap_new_row_clust.pdf)
  - **Figure 4**
    - [GO plots](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/WGCNA/GO_analysis/Parent_by_mod)
    - [Eigengene plots](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/WGCNA/eigengene_plots)

**glmmSeq differential expression and frontloading analysis**
- [glmmSeq script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/glmmSeq.Rmd)
- [Frontloading script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/Frontloading.Rmd)
- GO Enrichment Scripts
  - [DE by Origin](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/DEG_Enrich.Rmd)
  - [DE by Treatment](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/DEG_Trt_Enrich.Rmd)
  - [DE by Interaction](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/DEG_Int_Enrich.Rmd)
  - [Frontloaded](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/Frontloaded_Enrich.Rmd)
- [Outputs](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/glmmseq)
  - [Differential expression results, all genes](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/glmmseq/signif_genes.csv)
  - [Differentially expressed genes](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/glmmseq/DEGs)
  - [Frontloaded genes](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/glmmseq/frontloaded_genes.csv)
  - [GO Enrichment for DEGs](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/glmmseq/GOSeq)
    - [DE by Origin](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/glmmseq/GOSeq/GOseq_DEG_Origin.csv)
    - [DE by Treatment](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/glmmseq/GOSeq/Treatment/GOseq_DEG_Treatment.csv)
    - [DE by Interaction](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/glmmseq/GOSeq/Interaction/GOseq_DEG_Interaction.csv)
  - [GO Enrichment Frontloaded](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/glmmseq/Frontloaded_GOSeq)
    - [Results](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/glmmseq/Frontloaded_GOSeq/GOseq_Frontloaded.csv)

**Analysis of biomineralization-related genes**
- [BLAST Script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/glmmSeq/analysis/Biomineralization_Toolkit_Analysis.Rmd)
  - [BLAST Results](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/Biomineralization_blast_results.txt)
  - [Best hit table](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/Biomin_blast_Pocillopora_acuta_best_hit.csv)
- [Expression Plotting Script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/WGCNA/Biomineralization-toolkit-expression.Rmd)
- Outputs
  - [**Figure 5**](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/output/glmmseq/frontloaded_figures)
  - [**Figure 6**](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/output/WGCNA/biomin/biomin_Brown_Red_nonMod_compare_figs.pdf)

**Biomineralization-related genes ([Scucchia et al., 2021](https://doi.org/10.1111/gcb.15812); [Scucchia et al., 2021](https://doi.org/10.1098/rspb.2021.0328))**
- [Excel file](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/data/Biomineralization_Toolkit_FScucchia/Biomineralization_Toolkit_FScucchia.xlsx)
- [Fasta file](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/data/Biomineralization_Toolkit_FScucchia/Biomineralization_Toolkit_FScucchia.fasta)

**Reference Genome and annotations:**   
- Version 2 from: http://cyanophora.rutgers.edu/Pocillopora_acuta/

**NCBI Sequence uploads:**
- Bioproject: PRJNA934298  