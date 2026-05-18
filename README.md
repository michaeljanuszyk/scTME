# scTME
# This code was developed by the Longaker Laboratory at Stanford University (https://www.longakerlab.com)

Please contact Michael Januszyk directly with any questions or concerns

The following R and Python software versions used are listed below:
R v4.2.2
Seurat v4.9.9.9050
BPCells v0.1.0
SingleR v2.0.0
scATOMIC v2.0.2
zellkonverter v1.8.0
nichenetr v2.0.4
biomaRt v2.54.1
CellChat v1.5.0
infercnv v1.14.2
scMetabolism v0.2.1
SCENIC v1.3.1
TCGAbiolinks v3.18
TopHat v2.0.13
Python 3.9.19
scanpy v1.10.1
scrublet v0.2.3
scib v1.1.5

We recommend using a large size server to run this code. It took us roughly 6 weeks to run our code using all datasets (studyList.txt) on Dual AMD EPYC 7713 processors(2.0 GHz, 64-cores/128-threads) with 2TB of NVMe SSD memory. We believe on average the runtime scales in an O(nlogn) fasion, meaning that using only 50% of the samples would require considerably less than 50% of the runtime.

To run this code in cull, first download publicly available GEO or ArrayExpress datasets and convert each sample from each study into to a separate Seurat object. A list of the datasets used in the accompanying manuscript can be found in datasets_human.csv. Save each object as the accession code + .object.rds (e.g.,  GSM7114948.object.rds) and save them into folders named by study into a parent folder called readyForSeurat. For example, readyForSeurat/GSE190597/GSM7114948.object.rds. Next, please select either tmeMetaAnalysis_human.R or tmeMetaAnalysis_mouse.R, depending on the species. We recommend running this code in chunks, delineated by commented headings. The first few lines of this code will load tmeFunctions.R, which contains important functions used througout our code. As stated above, the expected run-time for this code is weeks to months depending on the hardware employed. 

Instructions for running this code on a very small number of samples can by found in exampleSetup.R, including expected output. Running the code in tmeMetaAnalysis_human.R on these samples takes approximately four hours.





