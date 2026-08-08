# Load required packages
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)

# Run SCENIC analysis on scRNA-seq Seurat object (objAll)
# 
# Note: Before running you must download the appropriate cisTarget databases (e.g., hg38 refseq-r80) 
# and place them in your relative db_dir.
objSce <- runScenicAnalysis(obj = objAll, org = "hgnc", db_dir = "data/cisTarget")

# Make SCENIC plots placing cell-level regulon data into original object (objAll)
# and produce the resulting FeaturePlots
makeScenicPlots( objSce, objAll, name )


