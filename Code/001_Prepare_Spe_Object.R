# The purpose of this script is to generate a spatialExperiment object (imcRtools),
# and cytoimageList objects (cytomapper) for images and masks to begin to facilitate
# the singleCell analysis.

# At this point, the input will be a Steinbock output folder.
# Images were preprocessed in two batches.
# The batches were manually combined into one steinbock folder (combined_steinbock)
# by manually appending images in the /img/ folder and likewise for /masks_deepcell/

################################################################################

# Clear environment
rm(list = ls())

# Load libraries
library(imcRtools) # For spatialExperiment object generation
library(cytomapper) # For plotting pixels and cells
library(tidyverse) # Data manipulation
library(RColorBrewer) # Color palettes
library(viridis) # Color palettes
library(scater) # Single-cell analysis
library(harmony) # Batch correction
library(BiocSingular) # Batch correction
library(patchwork) # Plotting
library(cowplot) # Plotting
library(dittoSeq) # Single-cell analysis
library(Rphenograph) # For clustering
library(igraph) # For clustering
library(ggplot2) # For plotting
library(ggpubr) # For statistical testing on plots

# Initialize a function to generate subdirectories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Prepare a subdirectory in the data directory to hold the pre-processed spe object
dataDir <- newDir("Data/001_Preprocessed")

# Prepare an output subdirectory
opDir <- newDir("Outputs/001_Prepare_Spe_Objects_Outputs")

# Prepare a subdirectory in the data directory to hold the pre-processed images
imgDir <- newDir("Data/Images")

#################################################################################
# Read in the Steinbock folder
spe <- read_steinbock("Data/combined_steinbock/")

# Make sample_id more interpretable
spe$sample_id <- sub("Slide1-6_CRCTMA1_CRCTMA2_", "ROI_", spe$sample_id)
spe$sample_id <- sub("Slide2-5_CRCTMA3_CRCTMA4_", "ROI_", spe$sample_id)

# Check that there are 26 unique sample_id
length(unique(spe$sample_id))

# For any markers that have a hyphen in their name, replace the hypen with an underscore
rownames(spe) <- gsub("\\-", "_", rownames(spe))

# Specify a subset of the channels that will be used for all downstream analysis
rowData(spe)$use_channel <- !grepl("DNA1|DNA2|ICSK1|ICSK2|ICSK3|CD14|CD31|CD11c|CD45RA|Perforin|Pan_actin|CD45RA", rownames(spe))

#################################################################################

# Metadata processing
metadata <- read.csv("Data/Updated_Metadata.csv")

# Add a column for Indication
spe$Indication <- metadata$Indication[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Patient ID
spe$Patient_ID <- metadata$Patient_ID[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Source
spe$Source <- metadata$Source[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Sex
spe$Sex <- metadata$Sex[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Age
spe$Age <- metadata$Age[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Diagnosis
spe$Diagnosis <- metadata$Diagnosis[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Stage
spe$Stage <- metadata$Stage[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Metastasis
spe$Mets <- metadata$Mets[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for "Keep"
spe$Keep <- metadata$Analysis[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for diagnosis
spe$Diagnosis <- metadata$Diagnosis[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for stage
spe$Stage <- metadata$Stage[match(spe$sample_id, metadata$Sample_ID)]

# Subset the spe object based on the "keep" vector
spe <- spe[,spe$Keep == "Keep"]
unique(spe$sample_id)

# Set cell names
colnames(spe) <- paste0('cell', seq_len(ncol(spe)))

# Filter the metadata
meta <- metadata[metadata$Analysis == "Keep",]
rownames(meta) <- meta$Sample_ID

################################################################################

# Initialize an empty list to hold the color vectors in the spe object
color_vectors <- list()

# Set colors for EOCRC and LOCRC
indication_colors <- c("EOCRC" = "firebrick2",
                       "LOCRC" = "cyan3")

sex_colors <- c("M" = "black",
                "F" = "grey")

source_colors <- c("JDH" = "green",
                   "WBH" = "purple")

# Append vectors to color list
color_vectors$Indication <- indication_colors
color_vectors$Sex <- sex_colors
color_vectors$Source <- source_colors

# Append color list to metadata slot in spe
metadata(spe)$color_vectors <- color_vectors

#################################################################################

# Read in the Images
img <- loadImages("Data/combined_steinbock/img/")

# Adjust sample IDs for images/masks to match spe
names(img) <- sub("Slide1-6_CRCTMA1_CRCTMA2_", "ROI_", names(img))
names(img) <- sub("Slide2-5_CRCTMA3_CRCTMA4_", "ROI_", names(img))

# Set a vector of the samples that we are keeping (only WBH and only tumor)
keep <- metadata[metadata$Analysis == "Keep",]$Sample_ID

# Filter the images based on the vector of samples to keep
img <- img[names(img) %in% keep]

# Read in the masks
masks <- loadImages("Data/combined_steinbock/masks_deepcell/", as.is = TRUE)

# Adjust sample IDs for images/masks to match spe
names(masks) <- sub("Slide1-6_CRCTMA1_CRCTMA2_", "ROI_", names(masks))
names(masks) <- sub("Slide2-5_CRCTMA3_CRCTMA4_", "ROI_", names(masks))

# Filter the masks based on the vector of samples to keep
masks <- masks[names(masks) %in% keep]

# Ensure that channel names for images match the spe object
channelNames(img) <- rownames(spe)

# Check that all images/masks are there and sample_ids are the same as spe
length(names(img))
names(img)
length(names(masks))
names(masks)

# Set image names across images/masks and spe
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img))

# Get metadata in the same order as the images
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[names(img),]

# Only keep the samples we want to include
metadata <- metadata[metadata$Analysis == "Keep",]
metadata <- metadata[names(img),]

# Generate vectors to append metadata to images
EorL <- metadata$Indication
Pat <- metadata$Patient_ID
Sex <- metadata$Sex
Diag <- metadata$Diagnosis
Stag <- metadata$Stage
Age <- metadata$Age
Source <- metadata$Source

# Assign the metadata to the images and masks as well
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img),
                                        Indication = EorL,
                                        Patient = Pat,
                                        Sex = Sex,
                                        Diagnosis = Diag,
                                        Stage = Stag, 
                                        Age = Age,
                                        Source = Source)

# Add a label slot to mcols for downstream plotting
mcols(img)$Label <- paste(paste(paste(mcols(img)$sample_id, mcols(img)$Sex, sep = " - "), mcols(img)$Source, sep = "\n"),
                          mcols(img)$Indication, sep = " - ")
mcols(masks)$Label <- paste(paste(paste(mcols(masks)$sample_id, mcols(masks)$Sex, sep = " - "), mcols(masks)$Source, sep = "\n"),
                          mcols(masks)$Indication, sep = " - ")

################################################################################

# Transform the counts to expression values
assay(spe, "exprs") <- asinh(counts(spe)/1)

# Run Dimensionality reduction on the spe object
set.seed(030161999)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs")

# Run PCA on the data to prepare for Harmony batch correction
set.seed(03061999)
spe <- runPCA(spe,
              subset_row = rownames(spe)[rowData(spe)$use_channel],
              exprs_values = "exprs",
              ncomponents = 50,
              BSPARAM = ExactParam())

# Run Harmony correction
set.seed(03061999)
out <- RunHarmony(spe, group.by.vars = c("sample_id", "Source"))
reducedDim(spe, "Harmony") <- reducedDim(out, "HARMONY")
spe <- runUMAP(spe, dimred = "Harmony", name = "UMAP_Harmony")

# Assess the batch correction
dittoDimPlot(spe, var = "Diagnosis", 
             reduction.use = "UMAP_Harmony", size = 0.2) +
  ggtitle("Patient ID on UMAP before correction")

################################################################################

# Set up a subdirectory to hold RPhenogrpah cluster results
PGdir <- newDir("Outputs/001_Prepare_Spe_Objects_Outputs/Phenograph_Clusters")

# Extract the Harmony corrected feature table from the spe object
mat <- reducedDim(spe, "Harmony")

# Set seed and run Rphenograph algorithm
set.seed(03061999)
out <- Rphenograph(mat, k = 45)

# Factor the clusters
clusters <- factor(membership(out[[2]]))

# Append clusters to the spe object
spe$pg <- clusters

# Remove the two samples we don't want from the anaysis
# Clustering still considered these samples since we clustered on the Harmony feature table which still
# includes these two samples. Therefore, removing them now shouldn't be a problem
spe <- spe[,spe$sample_id != "ROI_032"]
spe <- spe[,spe$sample_id != "ROI_035"]

# Take the image mean 
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$pg,
                                   statistics = "mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)

# Save image mean data as a csv file
imageMeanData <- as.data.frame(assay(image_mean, "exprs"))
# write.csv(imageMeanData, file = "Outputs/001_Prepare_Spe_Objects_Outputs/Phenograph_Clusters/PG_Image_Mean_Expression_Data.csv")

# PLot as a heatmap (max scaled)
PG_clusters <- dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", 
             cluster_cols = TRUE, 
             scaled.to.max = TRUE,
             heatmap.colors.max.scaled = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),  
             annot.by = c("pg"),
             show_colnames = TRUE)
ggsave("Outputs/001_Prepare_Spe_Objects_Outputs/Phenograph_Clusters/PG_Max_Scaled_Heatmap.tiff", PG_clusters, width = 8, height = 8, dpi = 300)

################################################################################

# Remove the two images we don't want from the analysis images and masks
img[["ROI_032"]] <- NULL
img[["ROI_035"]] <- NULL
masks[["ROI_032"]] <- NULL
masks[["ROI_035"]] <- NULL

# Save images and masks as Rds
saveRDS(img, file = "Data/Images/Final_Image_Dataset.Rds")
saveRDS(masks, file = "Data/Images/Final_Mask_Dataset.Rds")

# Normalize the images
img <- cytomapper::normalize(img, separateImages = TRUE)
images <- cytomapper::normalize(img, inputRange = c(0,0.2))

# Save normalized images as Rds
saveRDS(images, file = "Data/Images/Normalized_Final_Image_Dataset.Rds")

################################################################################

# Save the clustered spe object
saveRDS(spe, file = "Data/001_Preprocessed/PG_Clustered_Spe.Rds")

################################################################################

# This following code chunk is to assess differences in abundances of the identified clusters
# between early and late-onset CRC.

# Create a subdirectory in the output directory to store the cluster violins
clusterDir <- newDir("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins")

# Create a subdirectory in clusterDir to specify that these are by Indication
indicationDir <- newDir("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins/By_Indication")

# Differential marker expression
clusterNames <- unique(spe$pg)

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
# Iterate over the early-onset samples
for(i in EOsamples) {
  patient_cells <- sum(early$sample_id == i)
  
  # Iterate over each unique cluster
  for (j in clusterNames) {
    
    # Calculate the sum of each cluster for the sample
    cluster_cells <- sum(early$sample_id == i & early$pg == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append the proportion to the empty dataframe
    propDataE[i,j] <- proportion*100
  }
}

# Set column names of the proportion data to be the cluster names
colnames(propDataE) <- clusterNames

# Designate a new column labelled "Early"
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
# Iterate over the late-onset samples
for(i in LOsamples) {
  patient_cells <- sum(late$sample_id == i)
  
  # iterate over each unique cluster
  for (j in clusterNames) {
    
    # Calculate the sum of each cluster for the sample
    cluster_cells <- sum(late$sample_id == i & late$pg == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append the proportion to the empty dataframe
    propDataL[i,j] <- proportion*100
  }
}

# Set column names of the proportion data to be the cluster names
colnames(propDataL) <- clusterNames

# Designate a new column labelled as "Late"
propDataL$Group <- "Late"

# Combine the two results dataframes
propData <- rbind(propDataE, propDataL)

# Create a column for sample name
propData$Sample <- rownames(propData)

# Pivot longer the dataframe
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
  fileName <- paste("cluster", i, "By_Indication.tiff", sep = "_")
  ggsave(paste("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins/By_Indication/", fileName, sep = "/"), props, width = 10, height = 10, dpi = 300)
}

################################################################################

# Similar to the above, here, we will assess the differences between clusters based on sex
# instead of indication.

# Create a subdirectory in clusterDir to specify that these are by Indication
sexDir <- newDir("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins/By_Sex")

# Differential marker expression
clusterNames <- unique(spe$pg)

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
  
  # Iterate over the clusters
  for (j in clusterNames) {
    
    # Calculate the number of each cluster per sample
    cluster_cells <- sum(male$sample_id == i & male$pg == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append proportion to empty dataframe
    propDataM[i,j] <- proportion*100
  }
}

# Set column names to be cluster names
colnames(propDataM) <- clusterNames

# Add group column
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
  
  # Iterate over the cluster names
  for (j in clusterNames) {
    
    # Calculate the number of cells per cluster per sample
    cluster_cells <- sum(female$sample_id == i & female$pg == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append proportion information to empty dataframe
    propDataF[i,j] <- proportion*100
  }
}

# Set column names to be cluster names
colnames(propDataF) <- clusterNames

# Add group column
propDataF$Group <- "Female"

# Combine the two results dataframes
propData <- rbind(propDataM, propDataF)

# Add a column for sample name
propData$Sample <- rownames(propData)

# Pivot longer the data frame
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
  ggsave(paste("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins/By_Sex/", fileName, sep = "/"), props, width = 10, height = 10, dpi = 300)
}

################################################################################

# Following the same logic as above, we will now assess differences between the cluster
# based on both sex and indication.

# Create a subdirectory in clusterDir to specify that these are by Indication
sexIndDir <- newDir("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins/By_Sex_and_Indication")

# Differential marker expression
clusterNames <- unique(spe$pg)

# Now by Sex and Indication
maleE <- spe[,spe$Sex == "M" & spe$Indication == "EOCRC"]
maleL <- spe[,spe$Sex == "M" & spe$Indication == "LOCRC"]
femaleE <- spe[,spe$Sex == "F" & spe$Indication == "EOCRC"]
femaleL <- spe[,spe$Sex == "F" & spe$Indication == "LOCRC"]

# Subset the early spatial experiment object
MEsamples <- unique(maleE$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(MEsamples), ncol = length(clusterNames))
propDataME <- as.data.frame(propData)
rownames(propDataME) <- MEsamples
colnames(propDataME) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in MEsamples) {
  patient_cells <- sum(maleE$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(maleE$sample_id == i & maleE$pg == j)
    proportion <- cluster_cells / patient_cells 
    propDataME[i,j] <- proportion*100
  }
}

# Add colnames
colnames(propDataME) <- clusterNames

# Add group identifier
propDataME$Group <- "Male_Early"

#-----

# Subset the early spatial experiment object
MLsamples <- unique(maleL$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(MLsamples), ncol = length(clusterNames))
propDataML <- as.data.frame(propData)
rownames(propDataML) <- MLsamples
colnames(propDataML) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in MLsamples) {
  patient_cells <- sum(maleL$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(maleL$sample_id == i & maleL$pg == j)
    proportion <- cluster_cells / patient_cells 
    propDataML[i,j] <- proportion*100
  }
}

# Add colnames
colnames(propDataML) <- clusterNames

# Add group identifier
propDataML$Group <- "Male_Late"

#-----

# Subset the early spatial experiment object
FEsamples <- unique(femaleE$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(FEsamples), ncol = length(clusterNames))
propDataFE <- as.data.frame(propData)
rownames(propDataFE) <- FEsamples
colnames(propDataFE) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in FEsamples) {
  patient_cells <- sum(femaleE$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(femaleE$sample_id == i & femaleE$pg == j)
    proportion <- cluster_cells / patient_cells 
    propDataFE[i,j] <- proportion*100
  }
}

# Add colnames
colnames(propDataFE) <- clusterNames

# Add group identifier
propDataFE$Group <- "Female_Early"

#-----

# Subset the early spatial experiment object
FLsamples <- unique(femaleL$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(FLsamples), ncol = length(clusterNames))
propDataFL <- as.data.frame(propData)
rownames(propDataFL) <- FLsamples
colnames(propDataFL) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in FLsamples) {
  patient_cells <- sum(femaleL$sample_id == i)
  
  for (j in clusterNames) {
    cluster_cells <- sum(femaleL$sample_id == i & femaleL$pg == j)
    proportion <- cluster_cells / patient_cells 
    propDataFL[i,j] <- proportion*100
  }
}

# Add colnames
colnames(propDataFL) <- clusterNames

# Add group identifier
propDataFL$Group <- "Female_Late"

#-----
# Combine the two results dataframes
propData <- rbind(rbind(rbind(propDataME, propDataML), propDataFE), propDataFL)
propData$Sample <- rownames(propData)
propData <- propData %>%
  pivot_longer(-c("Sample", "Group"), names_to = "Cluster", values_to = "Proportion")

# Set custom colors
custom_colors <- c("Male_Early" = "red", "Male_Late" = "pink", "Female_Early" = "blue", "Female_Late" = "skyblue")

# Plot the proportions for each cluster
for (i in clusterNames) {
  subset <- propData[propData$Cluster == i,]
  
  runNo <- paste("Cluster", i, sep = " ")
  statement <- paste("Running preliminary stats on ", runNo, "\n", sep = "")
  cat(statement)
  
  # Assess equal variance
  cat("Levene test for equal variance \n")
  leveneRes <- leveneTest(Proportion ~ Group, data = subset)
  print(leveneRes)
  
  # Assess normality
  cat("SW test for normality \n")
  SW_res <- shapiro.test(subset$Proportion)
  print(SW_res)
  
  # Plot the violins with the statisical testing
  props <- ggplot(subset, aes(x = Group, y = Proportion, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(outliers = FALSE, width = 0.1) +
    geom_point() +
    theme_bw() +
    stat_pwc(method = "tukey_hsd") +
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
  ggsave(paste("Outputs/001_Prepare_Spe_Objects_Outputs/Cluster_Violins/By_Sex_and_Indication", fileName, sep = "/"), props, width = 10, height = 10, dpi = 300)
}

