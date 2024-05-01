library( dplyr); library(Seurat); library(ggplot2); library(Matrix); library(lsa)

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

 




# Log normalize
# Take top 2000 variable features
# Weight cosine distances by inverse rank
# Add rank as 3rd element in files

#### Arguments:
####   data = SeuratObject/SCE/matrix  ## Counts preferable
####   clusterLabels = NULL list of clusters (optional; else will do per cell)
####   level = 1  (or 2 or 3)
  scTME <- function( data, clusterLabels = NULL, level = 1 ) {
    data = log( data + 1e-10 )
    #data = log( data + 1 )
    clusterLabels = as.character(clusterLabels)    
    print( 'Stage 1' )
    #if( level == 2 ) {
      file_har = readRDS('scTME/centroid_har_fine_hvg.rds'); file_mnn = readRDS('scTME/centroid_mnn_fine_hvg.rds'); 
      file_scv = readRDS('scTME/centroid_scv_fine_hvg.rds')
    #}
    #mark_har = read.csv('scTME/harmony\ pbmc.markers.csv')
    #mark_mnn = read.csv('scTME/mnn\ pbmc.markers.csv')
    #mark_scv = read.csv('scTME/scvi\ pbmc.markers.csv')
#
    nm_har = file_har[[2]];    nm_mnn = file_mnn[[2]];    nm_scv = file_scv[[2]]

    centroid_har = file_har[[1]]
    centroid_mnn = file_mnn[[1]]
    centroid_scv = file_scv[[1]]

    centroids = list()
    nm_train = names( centroid_har[[1]] )
    nm_testt = rownames(data)
    nm_commo = intersect( nm_train, nm_testt )
    ii_train = match( nm_commo, nm_train )
    ii_testt = match( nm_commo, nm_testt )


    data = data[ ii_testt, ]
    data = t( t(data)/rowSums(t(data)) )
    data = as.matrix( as.data.frame ( data ) ) 

    for( i in 1:length(centroid_har) ) {; centroid_har[[i]] = centroid_har[[i]][ii_train]; }
    for( i in 1:length(centroid_mnn) ) {; centroid_mnn[[i]] = centroid_mnn[[i]][ii_train]; }
    for( i in 1:length(centroid_scv) ) {; centroid_scv[[i]] = centroid_scv[[i]][ii_train]; }

    for( i in 1:length(centroid_har) ) {; centroid_har[[i]] = centroid_har[[i]]/sum(as.double(centroid_har[[i]])); }
    for( i in 1:length(centroid_mnn) ) {; centroid_mnn[[i]] = centroid_mnn[[i]]/sum(as.double(centroid_mnn[[i]])); }
    for( i in 1:length(centroid_scv) ) {; centroid_scv[[i]] = centroid_scv[[i]]/sum(as.double(centroid_scv[[i]])); }

    print( 'Stage 2' )
    centroids[[level]] = list( centroid_har, centroid_mnn, centroid_scv )
    if( length(clusterLabels) == 0 ) {; 
      ## If doing on a per-cell basis
      clusterLabels = colnames(data);
      candidateCentroids = data
   } else {
     candidateCentroids = matrix( 0, nrow = nrow(data), ncol = length(unique(clusterLabels)) )
     for( i in 1:length(unique(clusterLabels)) ) {; print(i)
       if( sum(clusterLabels %in% sort(unique(clusterLabels))[i]) == 1 ) {
         candidateCentroids[ , i] =  data[,clusterLabels %in% sort(unique(clusterLabels))[i]]
       } else {
         candidateCentroids[ , i ] = rowSums(data[,clusterLabels %in% sort(unique(clusterLabels))[i]])/sum(clusterLabels %in% sort(unique(clusterLabels))[i])
       }
     } 
   }
   print( 'Stage 3' )
   distMat = array( 99, c(length(clusterLabels), length(centroids[[level]]), 20 ) )
   for( i in 1:length(unique(clusterLabels)) ) {; 
     for( j in 1:length( centroids[[level]] ) ) {
       for( k in 1:length( centroids[[level]][[j]] ) ) {
         #distMat[i,j,k] = 1 - cosine( candidateCentroids[,i], centroids[[level]][[j]][[k]][ii_train] )
         distMat[i,j,k] = 1 - cosine( candidateCentroids[,i], centroids[[level]][[j]][[k]] )
       }
     }
   }
   print( 'Stage 4' )
   win = matrix( 0 , nrow= length(unique(clusterLabels)), ncol = length(centroids[[level]]) )
   for( j in 1:length(centroids[[level]]) ) {
     for( i in 1:length(unique(clusterLabels)) ) {
       dists = distMat[ i, j, ]
       best = which( dists == min(dists) )
       if( length(best) > 0 ) {; best = best[1]; } 
       win[ i, j ] = best
     }
   }
   print( 'Stage 5' )
   nm_win = win
   nm_win[,1] = nm_har[ win[,1] ]
   nm_win[,2] = nm_mnn[ win[,2] ]
   nm_win[,3] = nm_scv[ win[,3] ]
   win = nm_win
   final = matrix( 0, nrow = length(unique(clusterLabels)) )
   for( i in 1:length(unique(clusterLabels)) ) {
     res = find_mode( win[i,] )
     if( length(res) > 1 ) {; res = res[1]; }
     final[i] = res
   }
   if( level==1 ) {
    final[ final=="Monocyte/macrophage"] <- "Myeloid"
    final[ final=="Plasma"] <- "Lymphoid"
    final[ final== "Mast"] <- "Myeloid"
    final[ final=="Pancreatic acinar"] <- "Epithelial"
    final[ final=="Hepatocyte"] <- "Epithelial"
    final[ final=="Retinal"] <- "Nerve"
    final[ final=="Alveolar"] <- "Epithelial"
    final[ final=="Neural crest"] <- "Epithelial"
    final[ final=="Keratinocyte"] <- "Epithelial"
    final[ final=="Hepatocyte"] <- "Epithelial"
    final[ final=="Renal tubule"] <- "Epithelial"
    final[ final=="Neutrophil"] <- "Myeloid"
    final[ final=="Dendritic"] <- "Myeloid"
    final[ final=="T cell"] <- "Lymphoid"
    final[ final=="B cell"] <- "Lymphoid"
  }
  stuff = clusterLabels
  for( i in 1:length(unique(clusterLabels)) ) {
    stuff[ clusterLabels %in% unique(clusterLabels)[i] ] = final[i]
  }
  return(stuff)
}



# temp4 = temp2[ as.integer( obj$seurat_clusters ) ]
