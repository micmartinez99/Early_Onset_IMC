# The purpose of this script is to visualize the clusters on the IMC images to aid
# With the phenotyping

# https://github.com/micmartinez99/Early_Onset_IMC

# To get this code onto your computer, simply clone the repository to your machine.
# You will have to go through and modify all the file paths. 

# Clear the environment
rm(list = ls())

# Load libraries
library(imcRtools)
library(cytomapper)
library(ggplot2)
library(viridis)
library(cowplot)
library(patchwork)

# Initialize a function to generate subdirectories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Create an output directory
opDir <- newDir("Outputs/002_Cluster_Visualization_Outputs")

################################################################################

# Read in the pre-processed spatialExperiment object
spe <- readRDS("Data/001_Preprocessed/PG_Clustered_Spe.Rds")

# Add Stage to spe object


sampleUMAP <- dittoDimPlot(spe, var = "Stage", 
                       reduction.use = "UMAP_Harmony", size = 0.2) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold"),
        legend.position = "bottom") +
  labs(title = "")
ggsave("Outputs/002_Cluster_Visualization_Outputs/Indication_Corrected_UMAP.tiff", sampleUMAP, width = 8, height = 8, dpi = 300)



# Visualize the clusters on the low dimensional embedding
pgUMAP <- dittoDimPlot(spe, var = "pg", 
             reduction.use = "UMAP_Harmony", size = 0.2, do.label = TRUE) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 30, face = "bold")) +
  labs(title = "PG Clusters")
ggsave("Outputs/002_Cluster_Visualization_Outputs/PG_clustered_UMAP.tiff", pgUMAP, width = 8, height = 8, dpi = 300)

# Visualize the marker expression on UMAP
markers <- rownames(spe)[rowData(spe)$use_channel]
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP_Harmony", 
                                assay = "exprs", size = 0.2, list.out = TRUE, axes.labels.show = FALSE, show.axes.numbers = FALSE, show.grid.lines = FALSE,
                                theme_classic())
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis(option = "A") +
                      theme(axis.title = element_text(face = "bold")))
markerPlots <- plot_grid(plotlist = plot_list) 
ggsave("Outputs/002_Cluster_Visualization_Outputs/New_Markers_on_Harmony_UMAP.tiff", markerPlots, width = 10, height = 10, dpi = 300)


# Read in the images and masks
img <- readRDS("Data/Images/Normalized_Final_Image_Dataset.Rds")
masks <- readRDS("Data/Images/Final_Mask_Dataset.Rds")

# Check that the names of img and masks match
all(names(img) == names(masks))

# Create subsets of the images and masks based on EO and LO for visualization
eocrc <- unique(spe[,spe$Indication == "EOCRC"]$sample_id)
locrc <- unique(spe[,spe$Indication == "LOCRC"]$sample_id)

# Subset the masks
eoMasks <- masks[names(masks) %in% eocrc]
loMasks <- masks[names(masks) %in% locrc]

# Subset the images
eoImg <- img[names(img) %in% eocrc]
loImg <- img[names(img) %in% locrc]

################################################################################

# Create an output subdirectory to hold the pixel images
pixelsDir <- newDir("Outputs/002_Cluster_Visualization_Outputs/Pixel_Images")

# Create an output subdirectory to hold the cell images
cellsDir <- newDir("Outputs/002_Cluster_Visualization_Outputs/Cell_Expr_Images")

# Subset spatialExperiment object to just cells in cluster 1
C1 <- spe[,spe$pg == "1"]

# Create a vector of datasets to run code in a loop
dataSets <- c("eoImg", "loImg")

# Iterate through early and late
for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  # Plot pixels
  cluster1 <- plotPixels(image = image,
                         mask = mask,
                         object = C1,
                         scale = FALSE,
                         cell_id = "ObjectNumber", img_id = "sample_id",
                         colour_by = c("DNA1", "CD3", "FoxP3", "CD4"),
                         outline_by = "pg",
                         bcg = list(DNA1 = c(0,1,1),
                                    CD3 = c(0,1,1),
                                    FoxP3 = c(0,2,1),
                                    CD4 = c(0,2,1)),
                         colour = list(pg = c("1" = "white")),
                         image_title = list(text = mcols(image)$Label,
                                            position = "top",
                                            colour = "white",
                                            margin = c(5,5),
                                            font = 4,
                                            cex = 1),
                         legend = list(colour_by.title.cex = 1,
                                       colour_by.labels.cex = 1,
                                       colour_by.title.font = 4,
                                       margin = 50,
                                       outline_by.title.cex = 1,
                                       outline_by.labels.cex = 1,
                                       outline_by.title.font = 4),
                         thick = TRUE,
                         scale = 10,
                         return_plot = TRUE)
  outputPlot <- ggdraw(cluster1$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "Cluster_1.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Pixel_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


# Iterate through early and late
for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  

# Examine the expression levels of FoxP3
C1_FoxP3 <- plotCells(mask,
          object = C1, 
          cell_id = "ObjectNumber", 
          img_id = "sample_id",
          colour_by = "FoxP3",
          exprs_values = "exprs",
          colour = list(FoxP3 = rev(turbo(100))),
          image_title = list(text = mcols(mask)$Label,
                             position = "top",
                             colour = "white",
                             margin = c(5,5),
                             font = 4,
                             cex = 1.3),
          legend = list(colour_by.title.cex = 1,
                        colour_by.labels.cex = 1,
                        colour_by.title.font = 4,
                        margin = 50),
          return_plot = TRUE)
outputPlot <- ggdraw(C1_FoxP3$plot, clip = "on")

# Generate fileName and save plot
fileName <- paste(file, "Cluster_1_FoxP3.tiff", sep = "_")
ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Cell_Expr_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

# Assess the different lymphocyte clusters
lymphocytes <- spe[,spe$pg == "1" | spe$pg == "13" | spe$pg == "15"]


for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  
  # Examine the expression levels of FoxP3
  C1_FoxP3 <- plotCells(mask,
                        object = lymphocytes, 
                        cell_id = "ObjectNumber", 
                        img_id = "sample_id",
                        colour_by = c("CD3"),
                        exprs_values = "exprs",
                        colour = list(CD3 = rev(turbo(100))),
                        image_title = list(text = mcols(mask)$Label,
                                           position = "top",
                                           colour = "white",
                                           margin = c(5,5),
                                           font = 4,
                                           cex = 1.3),
                        legend = list(colour_by.title.cex = 1,
                                      colour_by.labels.cex = 1,
                                      colour_by.title.font = 4,
                                      margin = 50),
                        return_plot = TRUE)
  outputPlot <- ggdraw(C1_FoxP3$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "Lymphocyte_C1_C13_C15_CD3.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Cell_Expr_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


###############################################################################

# Iterate through early and late
for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  # Plot pixels
  cluster1 <- plotPixels(image = image,
                         mask = mask,
                         object = lymphocytes,
                         scale = FALSE,
                         cell_id = "ObjectNumber", img_id = "sample_id",
                         colour_by = c("DNA1", "FoxP3", "CD4", "CD8a"),
                         outline_by = "pg",
                         bcg = list(DNA1 = c(0,1,1),
                                    FoxP3 = c(0,2,1),
                                    CD4 = c(0,2,1),
                                    CD8a = c(0,1,1)),
                         colour = list(pg = c("1" = "purple",
                                              "13" = "cyan1",
                                              "15" = "cyan1")),
                         image_title = list(text = mcols(image)$Label,
                                            position = "top",
                                            colour = "white",
                                            margin = c(5,5),
                                            font = 4,
                                            cex = 1),
                         legend = list(colour_by.title.cex = 1,
                                       colour_by.labels.cex = 1,
                                       colour_by.title.font = 4,
                                       margin = 50),
                         legend = list(outline_by.title.cex = 1,
                                       outline_by.labels.cex = 1,
                                       outline_by.title.font = 4),
                         thick = TRUE,
                         scale = 10,
                         return_plot = TRUE)
  outputPlot <- ggdraw(cluster1$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "Lymphocte_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Pixel_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

################################################################################
# Assess the different lymphocyte clusters
undefined <- spe[,spe$pg == "9" | spe$pg == "17" | spe$pg == "6"]


for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  # Plot pixels
  cluster1 <- plotPixels(image = image,
                         mask = mask,
                         object = undefined,
                         scale = FALSE,
                         cell_id = "ObjectNumber", img_id = "sample_id",
                         colour_by = c("DNA1", "Ecadherin", "Vimentin"),
                         outline_by = "pg",
                         bcg = list(DNA1 = c(0,1,1),
                                    Ecadherin = c(0,2,1),
                                    Vimentin = c(0,2,1)),
                         colour = list(pg = c("9" = "white",
                                              "17" = "burlywood2",
                                              "6" = "orange1")),
                         image_title = list(text = mcols(image)$Label,
                                            position = "top",
                                            colour = "white",
                                            margin = c(5,5),
                                            font = 4,
                                            cex = 1),
                         legend = list(colour_by.title.cex = 1,
                                       colour_by.labels.cex = 1,
                                       colour_by.title.font = 4,
                                       margin = 50),
                         legend = list(outline_by.title.cex = 1,
                                       outline_by.labels.cex = 1,
                                       outline_by.title.font = 4),
                         thick = TRUE,
                         scale = 10,
                         return_plot = TRUE)
  outputPlot <- ggdraw(cluster1$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "Undefined_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Pixel_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}



################################################################################

# Assess the IEL cluster
IELs <- spe[,spe$pg == "19"]

for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  # Plot pixels
  cluster1 <- plotPixels(image = image,
                         mask = mask,
                         object = IELs,
                         scale = FALSE,
                         cell_id = "ObjectNumber", img_id = "sample_id",
                         colour_by = c("Ecadherin", "CD3", "CD8a"),
                         outline_by = "pg",
                         bcg = list(
                                    Ecadherin = c(0,2,1),
                                    CD3 = c(0,2,1),
                                    CD8a = c(0,1,1)),
                         colour = list(pg = c("19" = "white")),
                         image_title = list(text = mcols(image)$Label,
                                            position = "top",
                                            colour = "white",
                                            margin = c(5,5),
                                            font = 4,
                                            cex = 1),
                         legend = list(colour_by.title.cex = 1,
                                       colour_by.labels.cex = 1,
                                       colour_by.title.font = 4,
                                       margin = 50),
                         legend = list(outline_by.title.cex = 1,
                                       outline_by.labels.cex = 1,
                                       outline_by.title.font = 4),
                         thick = TRUE,
                         scale = 10,
                         return_plot = TRUE)
  outputPlot <- ggdraw(cluster1$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "IEL_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Pixel_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}



for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  
  # Examine the expression levels of FoxP3
  C1_FoxP3 <- plotCells(mask,
                        object = IELs, 
                        cell_id = "ObjectNumber", 
                        img_id = "sample_id",
                        colour_by = c("CD8a"),
                        exprs_values = "exprs",
                        colour = list(CD8a = turbo(100)),
                        image_title = list(text = mcols(mask)$Label,
                                           position = "top",
                                           colour = "white",
                                           margin = c(5,5),
                                           font = 4,
                                           cex = 1.3),
                        legend = list(colour_by.title.cex = 1,
                                      colour_by.labels.cex = 1,
                                      colour_by.title.font = 4,
                                      margin = 50),
                        return_plot = TRUE)
  outputPlot <- ggdraw(C1_FoxP3$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "IELs_C19_CD8a.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Cell_Expr_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

























################################################################################

# Let's assign preliminary cell type labels and see if the trend holds
spe$cellTypes <- recode(spe$pg, 
                        "1" = "Activated Tregs",
                        "2" = "Cytotoxic Cells",
                        "3" = "Cytotoxic Cells",
                        "4" = "Epithelium",
                        "5" = "Stroma",
                        "6" = "Undefined",
                        "7" = "Stroma",
                        "8" = "Epithelium",
                        "9" = "Undefined",
                        "10" = "Stroma",
                        "11" = "Proliferating Epithelium",
                        "12" = "Macrophage",
                        "13" = "CD8a+ T-Lymphocytes",
                        "14" = "Epithelium",
                        "15" = "CD8a+ T-Lymphocytes",
                        "16" = "Proliferating Epithelium",
                        "17" = "Undefined",
                        "18" = "Epithelium",
                        "19" = "IELs")

spe$cellTypes <- ifelse(spe$pg == "1", "Activated Tregs", 
                        ifelse(spe$pg == "2", "Cytotoxic Cells",
                               ifelse(spe$pg == "3", "Cytotoxic Cells",
                                      ifelse(spe$pg == "4", "Epithelium",
                                             ifelse(spe$pg == "5", "Stroma",
                                                    ifelse(spe$pg == "6", "Undefined", 
                                                           ifelse(spe$pg == "7", "Stroma",
                                                                  ifelse(spe$pg == "8", "Epithelium",
                                                                         ifelse(spe$pg == "9", "Undefined",
                                                                                ifelse(spe$pg == "10", "Stroma",
                                                                                       ifelse(spe$pg == "11", "Proliferating Epithelium",
                                                                                              ifelse(spe$pg == "12", "Macrophage",
                                                                                                     ifelse(spe$pg == "13", "CD8a+ T-lymphocytes",
                                                                                                            ifelse(spe$pg == "14", "Epithelium",
                                                                                                                   ifelse(spe$pg == "15", "CD8a+ T-lymphocytes",
                                                                                                                          ifelse(spe$pg == "16", "Proliferating Epithelium",
                                                                                                                                 ifelse(spe$pg == "17", "Undefined",
                                                                                                                                        ifelse(spe$pg == "18", "Epithelium", "IELs"))))))))))))))))))

cellTypeColors <- setNames(c("#3F1B03", "#F4AD31", "#894F36", "#1C750C", "#EF8ECC", 
                           "#6471E2", "#4DB23B", "grey", "#BF0A3D"),
                           c("Stroma", "Epithelium", "Proliferating Epithelium", "Macrophage", "Cytotoxic Cells",
                             "CD8a+ T-lymphocytes", "IELs", "Undefined", "Activated Tregs"))
metadata(spe)$color_vectors$cellTypes <- cellTypeColors

# Factor the cellTypes
spe$cellTypes <- factor(spe$cellTypes, levels = c("Stroma", "Epithelium",
                                                  "Proliferating Epithelium", "Macrophage",
                                                  "Cytotoxic Cells", "CD8a+ T-lymphocytes", 
                                                  "IELs", "Activated Tregs", "Undefined"))


# Visualize the phenotyped clusters on the low dimensional embedding
phenoUMAP <- dittoDimPlot(spe, var = "cellTypes", 
                       reduction.use = "UMAP_Harmony", size = 0.2, do.label = TRUE,
                       color.panel = cellTypeColors) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold")) +
  labs(title = "")
ggsave("Outputs/002_Cluster_Visualization_Outputs/CellTypes_clustered_UMAP.tiff", phenoUMAP, width = 8, height = 8, dpi = 300)


# Append stage information to metadata
metadata <- read.csv("Data/Updated_Metadata.csv")
samples <- unique(spe$sample_id)
metadata <- metadata[metadata$Sample_ID %in% samples,]
write.csv(metadata, "Data/Final_Metadata.csv")


# Add a column for Stage
spe$Stage <- metadata$Stage[match(spe$sample_id, metadata$Sample_ID)]

# Save RDS
saveRDS(spe, "Data/CellTyped_pg_clusters.Rds")

# Take the image mean 
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$pg,
                                   statistics = "mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)

# Factor the cellTypes
image_mean$cellTypes <- factor(image_mean$cellTypes, levels = c("Stroma", "Epithelium",
                                                  "Proliferating Epithelium", "Macrophage",
                                                  "Cytotoxic Cells", "CD8a+ T-lymphocytes", 
                                                  "IELs", "Activated Tregs", "Undefined"))

# Plot as a heatmap
cellTypes_clusters <- dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
                            assay = "exprs", 
                            cluster_cols = TRUE, 
                            scaled.to.max = TRUE,
                            heatmap.colors.max.scaled = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),  
                            annot.by = c("cellTypes", "ncells"),
                            show_colnames = TRUE,
                            annotation_colors = list(cellTypes = metadata(spe)$color_vectors$cellTypes,
                                                     ncells = plasma(100)))
ggsave("Outputs/002_Cluster_Visualization_Outputs/Preliminary_CellType_Definition_Heatmap.tiff", cellTypes_clusters, width = 8, height = 8, dpi = 300)


image_mean <- aggregateAcrossFeatures(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = rownames(spe)[rowData(spe)$use_channel],
                                   use.assay.type = "exprs")
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)



for (i in dataSets) {
  # Match the masks to the appropriate dataset
  if (i == "eoImg") {
    image = eoImg
    mask = eoMasks
    file = "Early"
  } else {
    image = loImg
    mask = loMasks
    file = "Late"
  }
  
  
  # Examine the expression levels of FoxP3
  C1_FoxP3 <- plotCells(mask,
                        object = spe, 
                        cell_id = "ObjectNumber", 
                        img_id = "sample_id",
                        colour_by = "cellTypes",
                        exprs_values = "exprs",
                        image_title = list(text = mcols(mask)$Label,
                                           position = "top",
                                           colour = "white",
                                           margin = c(5,5),
                                           font = 4,
                                           cex = 1.3),
                        legend = list(colour_by.title.cex = 2,
                                      colour_by.labels.cex = 2,
                                      colour_by.title.font = 4,
                                      margin = 50),
                        return_plot = TRUE)
  outputPlot <- ggdraw(C1_FoxP3$plot, clip = "on")
  
  # Generate fileName and save plot
  fileName <- paste(file, "CellTypes.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/Cell_Expr_Images", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


################################################################################

# Create a subdirectory in the output directory to store the cluster violins
clusterDir <- newDir("Outputs/002_Cluster_Visualization_Outputs/CellType_Violins")

# Create a subdirectory in clusterDir to specify that these are by Indication
indicationDir <- newDir("Outputs/002_Cluster_Visualization_Outputs/CellType_Violins/By_Indication")

# Differential marker expression
clusterNames <- unique(spe$cellTypes)

# Split spe into early and late
early <- spe[,spe$Indication == "EOCRC"]
late <- spe[,spe$Indication == "LOCRC"]

# Subset the early spatial experiment object
EOsamples <- unique(early$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(EOsamples), ncol = length(clusterNames))
propDataE <- as.data.frame(propData)
rownames(propDataE) <- EOsamples
colnames(propDataE) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in EOsamples) {
  patient_cells <- sum(early$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(early$sample_id == i & early$cellTypes == j)
    proportion <- cluster_cells / patient_cells 
    propDataE[i,j] <- proportion*100
  }
}

colnames(propDataE) <- clusterNames

propDataE$Group <- "Early"

# Repeat this process for the late samples
LOsamples <- unique(late$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(LOsamples), ncol = 1)
propDataL <- as.data.frame(propData)
rownames(propDataL) <- LOsamples
colnames(propDataL) <- "Late"

# Calculate the proportion of each cell type per patient
# Initialize an empty matrix
propData <- matrix(0, nrow = length(LOsamples), ncol = length(clusterNames))
propDataL <- as.data.frame(propData)
rownames(propDataL) <- LOsamples
colnames(propDataL) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in LOsamples) {
  patient_cells <- sum(late$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(late$sample_id == i & late$cellTypes == j)
    proportion <- cluster_cells / patient_cells 
    propDataL[i,j] <- proportion*100
  }
}

colnames(propDataL) <- clusterNames
propDataL$Group <- "Late"

# Combine the two results dataframes
propData <- rbind(propDataE, propDataL)
propData$Sample <- rownames(propData)
write.csv(propData, file = "Outputs/002_Cluster_Visualization_Outputs/CellType_Proportion_Data.csv")
propData <- propData %>%
  pivot_longer(-c("Sample", "Group"), names_to = "Cluster", values_to = "Proportion")

# Set colors for early and late
custom_colors <- c("Early" = "#1B9E77", "Late" = "#7570B3")

# Plot the proportions for each cluster
for (i in clusterNames) {
  subset <- propData[propData$Cluster == i,]
  
  # Plot the violins with the statisical testing
  props <- ggplot(subset, aes(x = Group, y = Proportion, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(outliers = FALSE, width = 0.1) +
    geom_point() +
    theme_bw() +
    stat_compare_means(method = "t.test", method.args = list(var.equal = FALSE), hjust = -1, vjust = -1, size = 8) +
    labs(title = i,
         y = "Expression",
         fill = "") +
    scale_fill_manual(values = custom_colors) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 24),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 20),
          title = element_text(size = 18, face = "bold"))
  fileName <- paste("cluster", i, "By_Indication.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/CellType_Violins/By_Indication", fileName, sep = "/"), props, width = 10, height = 10, dpi = 300)
}

################################################################################

# Create a subdirectory in clusterDir to specify that these are by Indication
sexDir <- newDir("Outputs/002_Cluster_Visualization_Outputs/CellType_Violins/By_Sex")

# Differential marker expression
clusterNames <- unique(spe$cellTypes)

# Split spe into early and late
male <- spe[,spe$Sex == "M"]
female <- spe[,spe$Sex == "F"]

# Subset the early spatial experiment object
maleSamples <- unique(male$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(maleSamples), ncol = length(clusterNames))
propDataM <- as.data.frame(propData)
rownames(propDataM) <- maleSamples
colnames(propDataM) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in maleSamples) {
  patient_cells <- sum(male$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(male$sample_id == i & male$cellTypes == j)
    proportion <- cluster_cells / patient_cells 
    propDataM[i,j] <- proportion*100
  }
}

colnames(propDataM) <- clusterNames

propDataM$Group <- "Male"

# Repeat this process for the late samples
femalesamples <- unique(female$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(femalesamples), ncol = 1)
propDataF <- as.data.frame(propData)
rownames(propDataF) <- femalesamples
colnames(propDataF) <- "Female"

# Calculate the proportion of each cell type per patient
# Initialize an empty matrix
propData <- matrix(0, nrow = length(femalesamples), ncol = length(clusterNames))
propDataF <- as.data.frame(propData)
rownames(propDataF) <- femalesamples
colnames(propDataF) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in femalesamples) {
  patient_cells <- sum(female$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(female$sample_id == i & female$cellTypes == j)
    proportion <- cluster_cells / patient_cells 
    propDataF[i,j] <- proportion*100
  }
}

colnames(propDataF) <- clusterNames
propDataF$Group <- "Female"

# Combine the two results dataframes
propData <- rbind(propDataM, propDataF)
propData$Sample <- rownames(propData)
propData <- propData %>%
  pivot_longer(-c("Sample", "Group"), names_to = "Cluster", values_to = "Proportion")

# Set colors for early and late
custom_colors <- c("Male" = "blue", "Female" = "red")

# Plot the proportions for each cluster
for (i in clusterNames) {
  subset <- propData[propData$Cluster == i,]
  
  # Plot the violins with the statisical testing
  props <- ggplot(subset, aes(x = Group, y = Proportion, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(outliers = FALSE, width = 0.1) +
    geom_point() +
    theme_bw() +
    stat_compare_means(method = "t.test") +
    labs(title = paste("Cluster", i, sep = " "),
         y = "Expression") +
    scale_fill_manual(values = custom_colors) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 20),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          title = element_text(size = 18, face = "bold"))
  fileName <- paste("cluster", i, "By_Sex.tiff", sep = "_")
  ggsave(paste("Outputs/002_Cluster_Visualization_Outputs/CellType_Violins/By_Sex", fileName, sep = "/"), props, width = 10, height = 10, dpi = 300)
}




################################################################################

# Fisher's exact test
# Stratify the metadata into male and female
males <- metadata[metadata$Sex == "M",]
females <- metadata[metadata$Sex == "F",]

# Create contingency tables
mStage <- table(males$Indication, males$Stage)
FT_male <- fisher.test(mStage)
mosaicplot(mStage,
           main = "Males",
           color = TRUE)

# Create contingency tables
fStage <- table(females$Indication, females$Stage)
FT_female <- fisher.test(fStage)








