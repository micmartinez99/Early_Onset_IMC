# The purpose of this script is to continue to analyze the imaging mass cytometry data.
# In the previous script (001_Prepare_Spe_Object.R), we generated a spatialExperiment object,
# cytoImageList objects, performed dimensionality reduction, batch correction, and clustering.
# Now, in this script, we will phenotype the identified clusters and perform a single-cell analysis
# based on these celltype identifications.

################################################################################

# Clear the environment
rm(list = ls())

# Load libraries
library(imcRtools)
library(cytomapper)
library(RColorBrewer)
library(scater)
library(ggplot2)
library(ggpubr)
library(viridis)
library(cowplot)
library(patchwork)
library(corrplot)
library(Hmisc)
library(tidyverse)
library(dittoSeq)
library(ComplexHeatmap)
library(magick)
library(scales)
library(CATALYST)


# Initialize a function to generate subdirectories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Create an output directory
opDir <- newDir("Outputs/002_Analysis")

# Read in the metadata
metadata <- read.csv("Data/Updated_Metadata.csv")
samples <- unique(spe$sample_id)
metadata <- metadata[metadata$Sample_ID %in% samples,]
write.csv(metadata, "Data/Final_Metadata.csv")

###############################################################################

# Within the opDir, make a dimensionality reduction folder
dimRed <- newDir("Outputs/002_Analysis/UMAPs")

# Read in the pre-processed spatialExperiment object
spe <- readRDS("Data/001_Preprocessed/PG_Clustered_Spe.Rds")

# Add a column for Stage
spe$Stage <- metadata$Stage[match(spe$sample_id, metadata$Sample_ID)]

# Plot a UMAP for sample_id for Harmony corrected UMAP
sampleUMAP <- dittoDimPlot(spe, var = "sample_id", 
                           reduction.use = "UMAP_Harmony", size = 0.2) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold"),
        legend.position = "bottom") +
  labs(title = "")
ggsave("Outputs/002_Analysis/UMAPs/1_Sample_ID_Harmony_Corrected_UMAP.tiff", sampleUMAP, width = 8, height = 8, dpi = 300)

# Plot a UMAP for indication for Harmony corrected UMAP
indiUMAP <- dittoDimPlot(spe, var = "Indication", 
                           reduction.use = "UMAP_Harmony", size = 0.2, 
                           color.panel = metadata(spe)$color_vectors$Indication) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold"),
        legend.position = "bottom") +
  labs(title = "")
ggsave("Outputs/002_Analysis/UMAPs/2_Indication_Harmony_Corrected_UMAP.tiff", indiUMAP, width = 8, height = 8, dpi = 300)

# Plot a UMAP for source for Harmony corrected UMAP
sourceUMAP <- dittoDimPlot(spe, var = "Source", 
                         reduction.use = "UMAP_Harmony", size = 0.2, 
                         color.panel = metadata(spe)$color_vectors$Source) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold"),
        legend.position = "bottom") +
  labs(title = "")
ggsave("Outputs/002_Analysis/UMAPs/3_Source_Harmony_Corrected_UMAP.tiff", sourceUMAP, width = 8, height = 8, dpi = 300)

# Plot a UMAP for stage for Harmony corrected UMAP
stageUMAP <- dittoDimPlot(spe, var = "Stage", 
                           reduction.use = "UMAP_Harmony", size = 0.2) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold"),
        legend.position = "bottom") +
  labs(title = "")
ggsave("Outputs/002_Analysis/UMAPs/4_Stage_Harmony_Corrected_UMAP.tiff", stageUMAP, width = 8, height = 8, dpi = 300)

# Plot the phenogrpah clusters on UMAP
# Visualize the clusters on the low dimensional embedding
pgUMAP <- dittoDimPlot(spe, var = "pg", 
                       reduction.use = "UMAP_Harmony", size = 0.2, do.label = TRUE) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 30, face = "bold")) +
  labs(title = "PG Clusters")
ggsave("Outputs/002_Analysis/UMAPs/5_PG_Cluster_Nums_Harmony_Corrected_UMAP.tiff", pgUMAP, width = 8, height = 8, dpi = 300)

# Clean up the rownames
rownames(spe)[rowData(spe)$use_channel] <- c("aSMA", "Vimentin", "CD163", "Pan CK", "CD15",
                                             "CD45", "FoxP3", "CD4", "E-Cadherin", "CD68", "CD8a", 
                                             "GZMB", "Ki67", "Collagen I", "CD3", "CD45RO")

# Visualize the marker expression on UMAP
markers <- rownames(spe)[rowData(spe)$use_channel]
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP_Harmony", 
                                assay = "exprs", 
                                size = 0.2, 
                                list.out = TRUE, 
                                axes.labels.show = FALSE, 
                                show.axes.numbers = FALSE, 
                                show.grid.lines = FALSE,
                                theme_classic())

# Plot
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis(option = "B") +
                      theme(axis.title = element_text(face = "bold"),
                            axis.ticks = element_blank(),
                            strip.text = element_text(face = "bold"),
                            title = element_text(face = "bold")))
markerPlots <- plot_grid(plotlist = plot_list) 
ggsave("Outputs/002_Analysis/UMAPs/6_Markers_on_Harmony_Corrected_UMAP.tiff", markerPlots, width = 10, height = 10, dpi = 300)

###############################################################################

# Assign celltypes to each cluster based on UMAPs 
spe$cellTypes <- ifelse(spe$pg == "1", "Tregs", 
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
                                                                                              ifelse(spe$pg == "12", "M2 Macrophage",
                                                                                                     ifelse(spe$pg == "13", "T-Helpers",
                                                                                                            ifelse(spe$pg == "14", "Epithelium",
                                                                                                                   ifelse(spe$pg == "15", "CTLs",
                                                                                                                          ifelse(spe$pg == "16", "Proliferating Epithelium",
                                                                                                                                 ifelse(spe$pg == "17", "Undefined",
                                                                                                                                        ifelse(spe$pg == "18", "Epithelium", "IELs"))))))))))))))))))

# Add a color vector to the metadata slot for celltypes
cellTypeColors <- setNames(c("#3F1B03", "#894F20","#F4AD31", "#1C750C", "#EF8EC1", 
                             "#6471E1", "#F4800C", "cyan3", "darkgrey", "#BF0A3D"),
                           c("Stroma", "Epithelium", "Proliferating Epithelium", "M2 Macrophage", "Cytotoxic Cells",
                             "CTLs", "T-Helpers", "IELs", "Undefined", "Tregs"))
metadata(spe)$color_vectors$cellTypes <- cellTypeColors

# Factor the cellTypes
spe$cellTypes <- factor(spe$cellTypes, levels = c("Stroma", "Epithelium",
                                                  "Proliferating Epithelium", "M2 Macrophage",
                                                  "Cytotoxic Cells", "CTLs", "T-Helpers",
                                                  "IELs", "Tregs", "Undefined"))

# Visualize the phenotyped clusters on the low dimensional embedding
phenoUMAP <- dittoDimPlot(spe, var = "cellTypes", 
                          reduction.use = "UMAP_Harmony", 
                          size = 0.2, 
                          do.label = TRUE,
                          labels.size = 4,
                          labels.highlight = TRUE,
                          color.panel = cellTypeColors) +
  theme(axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 24, face = "bold")) +
  labs(title = "")
ggsave("Outputs/002_Analysis/UMAPs/7_CellTypes_Harmony_Corrected_UMAP.tiff", phenoUMAP, width = 8, height = 8, dpi = 300)

# Create a new colData entry in spe that will be a clean ID for plotting purposes
spe$ID <- paste(gsub("_", " ", spe$sample_id), spe$Indication, sep = " ")

# Read in the single cell masks 
masks <- readRDS("Data/Images/Final_Mask_Dataset.Rds")
mcols(masks)$ID <- paste(gsub("_", " ", mcols(masks)$sample_id), mcols(masks)$Indication, sep = " ")

# Pick an ROI to visualize (you can change this to any ROI in the dataset)
vis <- "ROI_002"
cur_masks <- masks[names(masks) %in% vis]

# Plot cell clusters on this ROI to visualize clustering
cellClusters <- plotCells(cur_masks,
          object = spe, 
          cell_id = "ObjectNumber", 
          img_id = "ID",
          colour_by = "cellTypes",
          colour = list(cellTypes = metadata(spe)$color_vectors$cellTypes),
          image_title = list(text = mcols(cur_masks)$ID,
                             position = "top",
                             colour = "white",
                             margin = c(5,5),
                             font = 2,
                             cex = 1.3),
          legend = list(colour_by.title.cex = 1,
                        colour_by.labels.cex = 1,
                        colour_by.title.font = 3,
                        margin = 50),
          return_plot = TRUE)
outputPlot <- ggdraw(cellClusters$plot, clip = "on")
ggsave("Outputs/002_Analysis/CellClusters_On_Tissues.tiff", outputPlot, width = 10, height = 10, dpi = 300)

################################################################################
# Heatmap analysis 

# Aggregate clusters across cells by taking the mean cluster counts per cell
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$pg,
                                   statistics = "mean",
                                   use.assay.type = "counts")

# Transform counts to expression values by arcsin transforming
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)

# Save the image mean data as a dataframe for manual plotting
imdata <- as.data.frame(t(assay(image_mean, "exprs")))

# Create a metadata dataframe for the PG clusters
clusterMeta <- data.frame(Cluster_Num = c(1:19),
                          Phenotype = c("Tregs", "Cytotoxic Cells",
                                        "Cytotoxic Cells", "Epithelium", "Stroma",
                                        "Undefined", "Stroma", "Epithelium", "Undefined",
                                        "Stroma", "Proliferating Epithelium", "M2 Macrophage",
                                        "T-Helpers", "Epithelium", "CTLs",
                                        "Proliferating Epithelium", "Undefined", "Epithelium", "IELs"))
clusterMeta$Cluster_Num <- paste("Cluster", rownames(clusterMeta), sep = " ")

# Factor the celltypes
clusterMeta$Phenotype <- factor(clusterMeta$Phenotype, levels = 
                                  c("Stroma", "Epithelium",
                                    "Proliferating Epithelium", "M2 Macrophage",
                                    "Cytotoxic Cells", "CTLs", "T-Helpers",
                                    "IELs", "Tregs", "Undefined"))


# Function to max scale the image mean data
max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply the function to the dataframe
maxScaled <- as.data.frame(lapply(imdata, max_scale))

# Clean up the rownames for the max-scaled data and for the metadata
rownames(maxScaled) <- paste("Cluster", rownames(maxScaled), sep = " ")
rownames(clusterMeta) <- as.character(rownames(maxScaled))
clusterMeta$Cluster_Num <- NULL

# Ensure the colors are in the order of the factored levels
phenotype_levels <- levels(clusterMeta$Phenotype)
phenotype_colors <- setNames(c("#3F1B03", "#F4AD31", "#894F20", "#1C750C", "#EF8EC1", "#6471E1", "#F4800C", "cyan3", "#BF0A3D", "grey"), 
                             phenotype_levels)

# Create the row annotation
CTAnno <- rowAnnotation(
  `Cell Type` = clusterMeta$Phenotype,
  col = list(`Cell Type` = phenotype_colors),
  show_annotation_name = FALSE
)

# Order the columns
colOrder <- c("Ecadherin", "Pan_keratin", "Ki67", 
              "aSMA", "Vimentin", "CollagenI",
              "CD45", "CD45RO", "CD3", "CD8a", "CD4", "FoxP3", 
              "CD68", "CD163", "CD15", "GranzymeB")
maxScaled <- maxScaled[,colOrder]

# Clean up column names
colnames(maxScaled) <- c("E-cadherin", "Pan-CK", "Ki67", "aSMA", "Vimentin", "Collagen I", "CD45",
                         "CD45RO", "CD3", "CD8a", "CD4", "FoxP3", "CD68", "CD163", "CD15", "GZMB")

# Set splitting pattern
hmSplit <- rep(1:5, c(3,3,6,2,2))

# Add metadata features
anno <- colData(image_mean) %>% 
  as.data.frame %>%
  select(cellTypes, ncells)

# Add number of cells per cluster as a row annotation
ha_meta <- rowAnnotation(`Number Cells` = anno_barplot(anno$ncells, width = unit(30, "mm")))

# Plot the heatmap                     
H <- Heatmap(as.matrix(maxScaled),
             show_column_names = TRUE,
             name = "Max Scale",
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             left_annotation = CTAnno,
             column_title = NULL,
             column_split = hmSplit,
             gap = unit(1, "mm"),
             row_title = NULL,
             col = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),
             row_dend_width = unit(2, "cm"),
             column_dend_height = unit(2, "cm"),
             row_dend_gp = gpar(lwd = 2.5),
             column_dend_gp = gpar(lwd = 2.5),
             right_annotation = ha_meta)
tiff("Outputs/002_Analysis/Publication_Heatmap.tiff", width = 10, height = 10, units = "in", res = 300)
draw(H)
dev.off()

# Alternative heatmap visualization
# Factor the cellTypes
image_mean$cellTypes <- factor(image_mean$cellTypes, levels = c("Stroma", "Epithelium",
                                                                "Proliferating Epithelium", "M2 Macrophage",
                                                                "Cytotoxic Cells", "CD8a+ T-lymphocytes", "CD4+ T-lymphocytes",
                                                                "IELs", "Tregs", "Undefined"))

# Aggregate clusters across cells by taking the mean cluster counts per cell
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$pg,
                                   statistics = "mean",
                                   use.assay.type = "counts")

# Transform counts to expression values by arcsin transforming
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)


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

################################################################################

# Create a new subdirectory for single cell analysis
scAnalysis <- newDir("Outputs/002_Analysis/SingleCell")

# Plot a compositional barplot
composition <- dittoBarPlot(spe, 
             var = "cellTypes", 
             group.by = "sample_id",
             split.by = "Indication",
             split.adjust = list(scales = "free")) +
  scale_fill_manual(values = metadata(spe)$color_vectors$cellTypes) +
  scale_y_continuous(labels = percent_format()) +
  labs(y = "Percent of Cells",
       x = "",
       fill = "",
       title = "") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 24),
        strip.text = element_text(face = "bold", size = 24),
        legend.position = "none")
ggsave("Outputs/002_Analysis/SingleCell/Compositional_Barplot_EOvsLO.tiff", composition, width = 8, height = 8, dpi = 300)
  
# Number of cellTyped cells per image
numCells <- dittoBarPlot(spe, 
             scale = "count",
             var = "cellTypes", 
             group.by = "sample_id",
             split.by = "Indication",
             split.adjust = list(scales = "free")) +
  scale_fill_manual(values = metadata(spe)$color_vectors$cellTypes) +
  labs(y = "Number of Cells",
       x = "",
       fill = "",
       title = "") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 24),
        strip.text = element_text(face = "bold", size = 24),
        axis.ticks.x = element_blank(),
        legend.position = "none")
ggsave("Outputs/002_Analysis/SingleCell/NumCells_Barplot_EOvsLO.tiff", numCells, width = 8, height = 8, dpi = 300)

# Save SPE in CATALYST-compatible object with renamed colData entries and 
# new metadata information
spe_cat <- spe 
spe_cat$sample_id <- factor(spe$sample_id)
spe_cat$condition <- factor(spe$Indication)
spe_cat$cluster_id <- factor(spe$cellTypes)

# Add celltype information to metadata
metadata(spe_cat)$cluster_codes <- data.frame(cellTypes = factor(spe_cat$cellTypes))

# Plot pseduobulk MDS plot
MDSplot <- pbMDS(spe_cat, 
      by = "cluster_id", 
      features = rownames(spe_cat)[rowData(spe_cat)$use_channel], 
      label_by = "cluster_id", 
      k = "cellTypes") +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$cellTypes) +
  theme(axis.title.x = element_text(face = "bold", size = 24),
        axis.title.y = element_text(face = "bold", size = 24),
        legend.position = "bottom") +
  guides(color = FALSE, 
         size = guide_legend()) +
  labs(y = "MDS 2",
       x = "MDS 1",
       color = "",
       size = "N Cells")
ggsave("Outputs/002_Analysis/SingleCell/Celltype_pbMDS.tiff", MDSplot, width = 8, height = 8, dpi = 300)

################################################################################

# Exactly as we did in the previous script, we will now look at the differential
# abundance of each cluster now that they are phenotyped.

# Create a subdirectory in the output directory to store the cluster violins
clusterDir <- newDir("Outputs/002_Analysis/Cluster_Violins")

# Create a subdirectory in clusterDir to specify that these are by Indication
indicationDir <- newDir("Outputs/002_Analysis/Cluster_Violins/Indication")

# Get a vector of unique cluster names
clusterNames <- unique(spe$cellTypes)

# Split spe into early and late
early <- spe[,spe$Indication == "EOCRC"]
late <- spe[,spe$Indication == "LOCRC"]

# Subset the early spatial experiment object
EOsamples <- unique(early$sample_id)

# Initialize an empty matrix to hold proportion data
propData <- matrix(0, nrow = length(EOsamples), ncol = length(clusterNames))
propDataE <- as.data.frame(propData)
rownames(propDataE) <- EOsamples
colnames(propDataE) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in EOsamples) {
  # Calculate the number of cells per patient
  patient_cells <- sum(early$sample_id == i)
  
  # Iterate over the cluster names
  for (j in clusterNames) {
    
    # Sum the number of each cell type per sample and calculate proportion
    cluster_cells <- sum(early$sample_id == i & early$cellTypes == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append the proportion to the empty dataframe
    propDataE[i,j] <- proportion*100
  }
}

# Set column names of proportion data to be cluster names
colnames(propDataE) <- clusterNames

# Set a new Group column as "EOCRC"
propDataE$Group <- "EOCRC"

# Repeat this process for the late samples
LOsamples <- unique(late$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(LOsamples), ncol = 1)
propDataL <- as.data.frame(propData)
rownames(propDataL) <- LOsamples
colnames(propDataL) <- "LOCRC"

# Calculate the proportion of each cell type per patient
# Initialize an empty matrix
propData <- matrix(0, nrow = length(LOsamples), ncol = length(clusterNames))
propDataL <- as.data.frame(propData)
rownames(propDataL) <- LOsamples
colnames(propDataL) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in LOsamples) {
  
  # Calculate the number of cells per sample
  patient_cells <- sum(late$sample_id == i)
  
  # Iterate over the cluster names
  for (j in clusterNames) {
    
    # Calculate the number of each cell type per sample and calculate proportion
    cluster_cells <- sum(late$sample_id == i & late$cellTypes == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append proportion data to the empty dataframe
    propDataL[i,j] <- proportion*100
  }
}

# Set column names of proportion data to be cluster names
colnames(propDataL) <- clusterNames

# Add a group column and set as "LOCRC"
propDataL$Group <- "LOCRC"

# Combine the two results dataframes
propData <- rbind(propDataE, propDataL)

# Add a column for sample ID
propData$Sample <- rownames(propData)

# Save results as a csv file
write.csv(propData, file = "Outputs/002_Analysis/Cluster_Violins/CellType_Proportion_Data.csv")

# Pivot longer the dataframe
propData <- propData %>%
  pivot_longer(-c("Sample", "Group"), names_to = "Cluster", values_to = "Proportion")

# Set colors for early and late
custom_colors <- c("EOCRC" = "firebrick2", "LOCRC" = "cyan3")
comparisons <- list(c("EOCRC", "LOCRC"))

# Plot the proportions for each cluster
for (i in clusterNames) {
  subset <- propData[propData$Cluster == i,]
  
  # Plot the violins with the statisical testing
  props <- ggplot(subset, aes(x = Group, y = Proportion, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(outliers = FALSE, width = 0.1) +
    geom_point() +
    theme_bw() +
    stat_compare_means(method = "t.test", comparisons = comparisons, label = "p.format",
                       method.args = list(var.equal = FALSE), vjust = -3, size = 8, 
                       bracket.size = 1,
                       bracket.nudge.y = 1) +
    labs(title = i,
         y = "Percent of Cell type / Image",
         fill = "") +
    scale_fill_manual(values = custom_colors) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 24),
          axis.text.y = element_text(size = 18),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 20),
          title = element_text(size = 18, face = "bold"))
  fileName <- paste("cluster", i, "By_Indication.tiff", sep = "_")
  ggsave(paste("Outputs/002_Analysis/Cluster_Violins/Indication", fileName, sep = "/"), props, width = 10, height = 10, dpi = 300)
}

################################################################################

# Create a data directory to hold the processed spe
dataDir <- newDir("Data/002_Processed")

# Save rds
saveRDS(spe, file = "Data/002_Processed/CellTyped_singleCell_Analysis_Spe.Rds")

################################################################################

