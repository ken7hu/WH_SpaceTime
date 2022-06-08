# WH_SpaceTime
custom code for analysis of data in Space-time Mapping Identifies Concerted Multicellular Patterns and Gene Programs in Healing Wounds and their Conservation in Cancers.
Each R script is reponsible for analysis and visualization of a subset of figure panels
######################################################################################

cytof_to_seurat.R : conversion of cytof data in tabular format to a Seurat object (Fig. 1E)

phemd_cytof.R : running phemd on the converted cytof Seurat object (Fig. 1E)

phemd_scrna_seq.R : phemd directly on scRNA-Seq Seurat object (Fig. 2J)

tileplot.R : generalized tileplot method for cell freqs, average gene expression within a population and average factor level used for (Fig. 2H,2I, 3C, 4C, 5A, 5B)

generalized_nmf.R: generalized function for computing NMF decomp on a given Seurat object. Used for (Fig. 5)

nmf_downstream_analysis.R: graphing utilities for visualizing average factor expression over space time and gene weight contributions for computed factors (Fig. 5A, 5B, 5E-5G, S5A,B,C)

correlation_matrices.circles.R : graphing utilities for generating correlograms and radial factor level movie (Movie S5) (Fig. 5C, S5F-G, S5H)

imaris_analysis.R :scrips for statistical analysis and plotting of imaris output. Used for (Fig. 7 and S7)
