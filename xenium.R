library(Seurat)
library(terra)

path <- "C:/Users/AlainaKK/Downloads/Xenium/output-XETG00074__0010499__IC071922__20231001__221931"
# Load the Xenium data
xenium_obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium_obj <- subset(xenium_obj, subset = nCount_Xenium > 0)

# Extract Cell Centroids and Molecule Coordinates
# The counts matrix (“matrix”): This contains expression data for cells and features.
# Cell centroids in pixel coordinate space (“centroids”): Provides cell centroid coordinates (x, y).
# Molecule pixel coordinates (“microns”): Gives the pixel coordinates of individual molecules.
xenium_read <- ReadXenium(path, outs = c("matrix", "microns"), type = "centroids", mols.qv.threshold = 20)


# Function to compute each cell's area and correlate it with clusters
get_cell_area_and_correlate <- function(xenium_obj) {
  
  xenium_obj[['area']] <- paste0("(",xenium_read[['centroids']][['x']],",",xenium_read[['centroids']][['y']],")")
  
  # Run UMAP to reduce dimensions for visualization
  xenium_obj <- SCTransform(xenium_obj, assay = "Xenium")
  xenium_obj <- RunPCA(xenium_obj, npcs = 30, features = rownames(xenium_obj))
  xenium_obj <- RunUMAP(xenium_obj, dims = 1:30)
  
  # Identify clusters of cells
  xenium_obj <- FindNeighbors(xenium_obj, reduction = "pca", dims = 1:30)
  xenium_obj <- FindClusters(xenium_obj, resolution = 0.3)
  
  # Plot UMAP, colored by cell area ???
  # Scatter plot
  FeaturePlot(xenium_obj, features = rownames(xenium_obj))
  # Violin plot
  VlnPlot(xenium_obj, features = rownames(xenium_obj))
  
  # Plot UMAP, colored by clusters
  DimPlot(xenium_obj)
  
  # correlation ???
  ClusterAreaCorr <- cor.test(xenium_obj[['area']], xenium_obj[['seurat_clusters']], method = "spearman")
  
  # Return the modified Seurat object, which now includes cell area metadata
  # and clustering information
  return(xenium_obj)
}

xenium_seurat <- get_cell_area_and_correlate(xenium_obj)








