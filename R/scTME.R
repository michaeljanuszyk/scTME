function (obj, clusters, ref, level = 1, max.cells = 50000, return.details = F) {
    options(Seurat.object.assay.version = "v3")
    require(SingleR); library(SingleR)
    require(SingleCellExperiment); library(SingleCellExperiment)
    if (class(obj) != "SingleCellExperiment") {
        sce = as.SingleCellExperiment(obj)
    }
    else {
        sce = obj
    }
    if (dim(sce)[2] > max.cells) {
        x = sce
        keep <- sample(Cells(x), max.cells)
        x$keep <- "false"
        x$keep[keep] <- "true"
        x.down <- subset(x, subset = keep == "true")
        sce = x.down
    }
    common <- intersect(rownames(sce), rownames(ref))
    sce <- sce[common, ]
    ref <- ref[common, ]
    if (level == 1) {
        labels = ref$integrated.annotation.abbr
    }
    else {
        labels = ref$integrated.annotation
    }
    pred <- SingleR(test = sce, ref = ref, labels = labels, de.method = "wilcox", 
        clusters = clusters)
    newLabels = pred$labels
    names(newLabels) <- levels(as.factor(clusters))
    tempLabels = "tbd"
    for (c in unique(clusters)) {
        tempLabels[clusters == c] = newLabels[c]
    }
    if (return.details) {
        returnList = list(tempLabels, pred)
        names(returnList) = c("cell.level.predictions", "cluster.level.details")
        return(returnList)
    }
    else {
        return(tempLabels)
    }
}

