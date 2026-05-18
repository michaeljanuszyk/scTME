# Example code for running on two datasets

# First download raw data from:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=267718
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=284205

# Next load these functions

# Convert .h5 files to Seurat objects using standardized naming
curate_h5 <- function( geoName ) {
  fileList <- list.files( paste0( 'rawData/', geoName ), '*h5' )
  for (fileName in fileList) {; print( fileName )
    pathName <- paste0( 'rawData/', geoName, '/', fileName )
    objName <- substring(fileName,0,10)
    mat <- Read10X_h5( pathName, use.names = T, unique.features = T )
    obj <- CreateSeuratObject( mat, min.cells = 3, min.features = 200  )
    obj$orig.ident <- objName
    obj$study <- geoName
    dir.create( paste0( 'readyForSeurat/', geoName, '/' ), suppressWarnings(F) )
    saveName <- paste0( 'readyForSeurat/', geoName, '/',objName,'.object.rds' )
    saveRDS( obj, saveName, compress = F )
  }
}

# Convert .mtx file sets to Seurat objects using standardized naming
curate_mtx <- function( geoName ) {
  fileList <- list.files(paste0("rawData/",geoName),"*matrix.mtx")
  for( fileName in fileList ) {; print(fileName)
    GSM_name <- substring(fileName,0,10)
    dir_name <- paste0("rawData/",geoName,"/")
    mat <- ReadMtx(
      mtx = paste0(dir_name, fileName),
      features = paste0(dir_name, substring(fileName,0,nchar(fileName)-10),"features.tsv"),
      cells = paste0(dir_name, substring(fileName,0,nchar(fileName)-10),"barcodes.tsv")
    )
    obj <- CreateSeuratObject( mat, min.cells = 3, min.features = 200  )
    obj$orig.ident <- GSM_name
    obj$study <- geoName
    dir.create(paste0("readyForSeurat/",geoName,"/"),suppressWarnings(F))
    saveName <- paste0("readyForSeurat/",geoName,"/",GSM_name,".object.rds")
    saveRDS(obj,saveName, compress = F)
  }
}



# Now run these functions to load specific GEO datasets
curate_mtx('GSE267718')
curate_h5( 'GSE284205')


# Next, run the code in tmeAnalysis_human.R on samples from these two datasets





