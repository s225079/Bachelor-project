library(dplyr)
library(tidyr)
library(pheatmap)

# Example data preparation
path_heat <- paste('C:\\Users\\Miras\\OneDrive - Danmarks Tekniske Universitet\\Skrivebord\\6. semester\\Bachelor\\Genome analysis\\Collective heatmap', sep = '')

df1 <- read.table(paste0(path_heat, "/A.r_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df2 <- read.table(paste0(path_heat, "/R.sal_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df3 <- read.table(paste0(path_heat, "/B_int_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df4 <- read.table(paste0(path_heat, "/S_copri_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df5 <- read.table(paste0(path_heat, "/B_ova_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df6 <- read.table(paste0(path_heat, "/B_long_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df7 <- read.table(paste0(path_heat, "/P_clar_family_substrate_counts.txt"), sep = "\t", header = TRUE)
df8 <- read.table(paste0(path_heat, "/I_but_family_substrate_counts.txt"), sep = "\t", header = TRUE)

# list of data frames with names for each strain
strain_data <- list(
  'B. intestini' = df3,
  'S. salivae' = df2,
  'A. rectalis' = df1,
  'S. copri' = df4,
  'B. ovatus' = df5,
  'B. infantis' = df6,
  'P. clara' = df7,
  'I. butyriciproducens' = df8
)

substrate_colors <- c(
  'alpha-glucan' = 'lightpink',
  'arabinan'='orchid2',
  'beta-galactan'='olivedrab3',
  'beta-glucan' = 'mediumpurple1',
  'beta-mannan'='royalblue',
  'fructan'='orange',
  'pectin' = 'lightcoral',
  'starch' = 'turquoise2',
  'xylan' = 'seagreen3')

# Combine all strains into a single data frame
combined_data <- bind_rows(strain_data, .id = "Strain")

# Create unique row identifiers for families and substrates
combined_data <- combined_data %>%
  mutate(RowID = paste(family, substrate, sep = "_")) 

# Remove missing values in substrates
combined_data <- combined_data %>% drop_na(substrate)

# Order by substrate
combined_data <- combined_data %>% arrange(substrate)

# Create heatmap data
heatmap_data <- combined_data %>%
  select(RowID, Strain, count) %>%
  pivot_wider(names_from = Strain, values_from = count, values_fill = list(count = 0))

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data %>% select(-RowID)) 
rownames(heatmap_matrix) <- heatmap_data$RowID  


strain_order <- c("A. rectalis", "S. salivae", "B. intestini", "S. copri", "B. ovatus", "B. infantis", "P. clara", "I. butyriciproducens")
heatmap_matrix <- heatmap_matrix[, strain_order, drop = FALSE]

# Transpose matrix to flip rows and columns
heatmap_matrix <- t(heatmap_matrix)

# Log normalization (log1p)
heatmap_matrix_normalized <- log1p(heatmap_matrix)


# Extract substrate names from RowID
substrate_names <- sub("^[^_]+_", "", colnames(heatmap_matrix))

# Create annotation data as a data frame
annotation_row <- data.frame(Substrate = substrate_names)

# Set row names to match heatmap row names
rownames(annotation_row) <- colnames(heatmap_matrix)

# Extract family names
family_names <- sub("_.*$", "", colnames(heatmap_matrix))


strain_labels <- c("italic('A. rectalis')", "italic('S. salivae')", "italic('B. intestini')",
                   "italic('S. copri')", "italic('B. ovatus')", "italic('B. infantis')",
                   "italic('P. clara')", "italic('I. butyriciproducens')")


# Convert strings to expressions
strain_labels_expr <- parse(text = strain_labels)

# heatmap 
pheatmap(heatmap_matrix_normalized, cluster_rows = TRUE, cluster_cols = FALSE,
         main = "CAZyme families linked to dietary fibers identified in bacterial strains",
         angle_col = 90, fontsize=15,
         annotation_col = annotation_row,
         annotation_colors = list(Substrate = substrate_colors),
         breaks = seq(min(heatmap_matrix_normalized), max(heatmap_matrix_normalized), length.out = 100),
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         border_color = NA,
         labels_col = family_names,
         labels_row = strain_labels_expr,
         legend = TRUE
)


