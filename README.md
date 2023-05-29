# WH_SpaceTime
custom code for analysis of data in Space-time Mapping Identifies Concerted Multicellular Patterns and Gene Programs in Healing Wounds and their Conservation in Cancers.
Each R script is reponsible for analysis and visualization of a subset of figure panels
######################################################################################
Merging_Cleaning_CD45pos.R : Merging of 5 libraries, processing and cleaning and subsetting of CD45 pos for subsequent drill down and reclustering to arrive at final CD45pos object in Figure 1

CD45neg_processing_manuscript : Takes the cd45neg subset of the merged cleaned dataset, reclusters and reUMAPs to generate figure 3A and the dotplot and tileplots in Supplementary Figure 3A, 3B

fibroblasts_processing_manuscript : Subsets out fibroblasts for Figure 3B-D and Supplementary Figures 3C-G including comparison to the differential genes sets from Buechler et al. 2021 (Figure S3E) and Vu et al. 2022 (Figure S3F).

cytof_to_seurat_phemd.R : conversion of cytof data in tabular format to a Seurat object followed by running phemd on the converted cytof Seurat object (Mendeley Extended Figure 1)

phemd_scrna_seq.R : phemd directly on scRNA-Seq CD45pos Seurat object (Fig. 1J)

tileplot.R : generalized tileplot method for cell freqs, average gene expression within a population and average factor level used for (1H, 1I, 2C, 3D, 3F, 3I, 4A, 4B)

monocle_analysis_mono_mac : analysis script for subsetting monomac mhcii low and hi populations for monocle analysis. Includes plotting of gene expression by pseudotime and trajectory superimposed on the UMAP dimensional reduction (Fig. 1D-G)

generalized_nmf.R: generalized function for computing NMF decomp on a given Seurat object. Used for (Fig. 4,5,6)

cellchat.R : script for running cellchat from Seurat objects and graphing utilities (Fig. S4D)

nmf_downstream_analysis.R: graphing utilities for visualizing average factor expression over space time and gene weight contributions for computed factors (Fig. 4A, 4B, 4E, 4F, 4G, and S4A, S4B)

correlation_matrices.circles.R : graphing utilities for generating correlograms and radial factor level movie (Fig. 4C, S4J)

nmf_backend_tumor_to_wh.R : script for comparing factors between wound healing and tumor contextx (Fig. 5,6 and S5,S6)

imaris_analysis.R :scrips for statistical analysis and plotting of imaris output. Used for (Fig. 3H, S3J, S4G, S4H,7 and S7)
