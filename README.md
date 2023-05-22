# WH_SpaceTime
custom code for analysis of data in Space-time Mapping Identifies Concerted Multicellular Patterns and Gene Programs in Healing Wounds and their Conservation in Cancers.
Each R script is reponsible for analysis and visualization of a subset of figure panels
######################################################################################
Merging_Cleaning_CD45pos.R : Merging of 5 libraries, processing and cleaning and subsetting of CD45 pos for subsequent drill down and reclustering to arrive at final CD45pos object in Figure 1

cytof_to_seurat_phemd.R : conversion of cytof data in tabular format to a Seurat object followed by running phemd on the converted cytof Seurat object  (Fig. 1E)

phemd_scrna_seq.R : phemd directly on scRNA-Seq Seurat object (Fig. 2J)

tileplot.R : generalized tileplot method for cell freqs, average gene expression within a population and average factor level used for (Fig. 2H,2I, 3C, 4C, 5A, 5B)

monocle_analysis_mono_mac : analysis script for subsetting monomac mhcii low and hi populations for monocle analysis. Includes plotting of gene expression by pseudotime and trajectory superimposed on the UMAP dimensional reduction (Fig. 2D-G)

generalized_nmf.R: generalized function for computing NMF decomp on a given Seurat object. Used for (Fig. 5)

cellchat.R : script for running cellchat from Seurat objects and graphing utilities (Fig. S5D and S5E)

nmf_downstream_analysis.R: graphing utilities for visualizing average factor expression over space time and gene weight contributions for computed factors (Fig. 5A, 5B, 5E-5G, S5A,B,C)

correlation_matrices.circles.R : graphing utilities for generating correlograms and radial factor level movie (Movie S5) (Fig. 5C, S5F-G, S5H)

nmf_backend_tumor_to_wh.R : script for comparing factors between wound healing and tumor contextx (Fig. 6 and S6)

imaris_analysis.R :scrips for statistical analysis and plotting of imaris output. Used for (Fig. 7 and S7)
