########## Initialization########
library(dplyr)
library(DT)
library(BiocManager)
library(shiny)
library(matrixStats)
library(tools)
library(curatedMetagenomicData)

library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(pheatmap)
library(stringr)
library(phyloseq)
library(SummarizedExperiment)
library(mia)
library(scater)
library(scales)
library(car)
library(ggpubr)
library(rstatix)

library(writexl)

########### Defining species##########

species_names <- c(
  '[Eubacterium] rectale' = 'A. rectalis',
  'Prevotella salivae'='S. salivae',
  'Prevotella copri'='S. copri',
  'Bacteroides ovatus' = 'B. ovatus',
  'Veillonella dispar'='V. dispar',
  'Bifidobacterium longum'='B. longum',
  'Paraprevotella clara' = 'P. clara',
  'Intestinimonas butyriciproducens' = 'I. butyriciproducens'
  
)



##### Collection base data #######

data <- sampleMetadata

dat_base <-sampleMetadata |>
  filter(age >= 18) |>
  filter(body_site == "stool") |>
  filter(disease == "healthy") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")

data_base <- colData(dat_base) %>% data.frame()

tax_base <- rowData(dat_base) %>% data.frame()

abu_base <- assay(dat_base) %>% data.frame()


# Filtration of abudance data

# Define the vector of names
names <- c("[Eubacterium] rectale", "Prevotella salivae", "Prevotella copri", 
           "Bacteroides ovatus", "Veillonella dispar", "Bifidobacterium longum", 
           "Paraprevotella clara", "Intestinimonas butyriciproducens")

# Filter the dataframe based on row names
abu_base_filtered <- abu_base[rownames(abu_base) %in% names, ]


n_base <- (ncol(abu_base_filtered)/nrow(data))*100

# Convert to long format for ggplot
abu_base_long <- abu_base_filtered %>%
  rownames_to_column(var = "Species") %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")



#Prevalence calculations

#Counting the number of columns=0
zero_counts <- numeric(nrow(abu_base_filtered))  

for (i in 1:nrow(abu_base_filtered)){
  zero_counts[i] <- sum(abu_base_filtered[i, ] == 0)
}


#Calculating the prevelance
num_columns <- ncol(abu_base_filtered)  

Prev <- numeric(nrow(abu_base_filtered))

for (i in 1:8){
  
  Prev[i] <- 100-((zero_counts[i]/num_columns)*100)
}

df_prev <- data.frame(Species=rownames(abu_base_filtered),Prevalence=Prev)
df_prev <- df_prev %>% mutate(category = "Healthy adults (all)")

########## Diets ##########

#Abundance omnivore diet
WesternStudy_n <-
  filter(sampleMetadata, age >= 18) |>
  filter(!is.na(diet)) |>
  filter(diet=="omnivore") |>
  filter(disease == "healthy") |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")

data_omnivore <- colData(WesternStudy_n) %>% data.frame()

n_omni <- (nrow(data_omnivore)/nrow(data))*100

abu_western_n <- assay(WesternStudy_n) %>% data.frame()
abu_western_n_filtered <- abu_western_n[rownames(abu_western_n) %in% names, ]


abu_western_n_long <- abu_western_n_filtered %>%
  rownames_to_column(var = "Species") %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")


#Prevalence omnivore

zero_counts_n <- numeric(nrow(abu_western_n_filtered)) 

for (i in 1:nrow(abu_western_n_filtered)){
  zero_counts_n[i] <- sum(abu_western_n_filtered[i, ] == 0)
}


#Calculating the prevelance
num_columns <- ncol(abu_western_n_filtered) 

Prev_west_n <- numeric(nrow(abu_western_n_filtered))


for (i in 1:1:nrow(abu_western_n_filtered)){
  
  Prev_west_n[i] <- 100-((zero_counts_n[i]/num_columns)*100)
}

df_west_n <- data.frame (Species=rownames(abu_western_n_filtered), Prevalence=Prev_west_n, category="Healthy adults (omnivore diet)")


#Abundance vegetarian
WesternStudy_y <-
  filter(sampleMetadata, age >= 18) |>
  filter(!is.na(diet)) |>
  filter(diet == "vegetarian") |>
  filter(disease == "healthy") |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")

data_vege <- colData(WesternStudy_y) %>% data.frame()

n_vege <- (nrow(data_vege)/nrow(data))*100

abu_western_y <- assay(WesternStudy_y) %>% data.frame()

abu_western_y_filtered <- abu_western_y[rownames(abu_western_y) %in% names, ]


abu_western_y_long <- abu_western_y_filtered %>%
  rownames_to_column(var = "Species") %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")


#Prevalence vegetarian
zero_counts_y <- numeric(nrow(abu_western_y_filtered))  

for (i in 1:nrow(abu_western_y_filtered)){
  zero_counts_y[i] <- sum(abu_western_y_filtered[i, ] == 0)
}


#Calculating the prevelance
num_columns <- ncol(abu_western_y_filtered) 

Prev_west_y <- numeric(nrow(abu_western_y_filtered))


for (i in 1:1:nrow(abu_western_y_filtered)){
  
  Prev_west_y[i] <- 100-((zero_counts_y[i]/num_columns)*100)}

df_west_y <- data.frame (Species=rownames(abu_western_y_filtered), Prevalence=Prev_west_y, category="Healthy adults (vegetarian diet)")

############ DIet #############



#collect data
IBD <-
  filter(sampleMetadata, age >= 18) |>
  filter(disease == "healthy") |>
  filter(body_site == "stool") |>
  filter(!is.na(diet)) |>
  filter(diet == "vegan") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")

data_IBD <- colData(IBD) %>% data.frame()


abu_IBD <- assay(IBD) %>% data.frame()


n_vegan <- (nrow(data_IBD)/nrow(data))*100


#Abundance
# Filter the dataframe based on row names
abu_IBD_filtered <- abu_IBD[rownames(abu_IBD) %in% names, ]


# Convert to long format 
abu_IBD_long <- abu_IBD_filtered %>%
  rownames_to_column(var = "Species") %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")


#Prevalence calculations (vegan)

#Counting the number of columns=0
zero_counts <- numeric(nrow(abu_IBD_filtered)) 


for (i in 1:nrow(abu_IBD_filtered)){
  zero_counts[i] <- sum(abu_IBD_filtered[i, ] == 0)
}
zero_counts



#Calculating the prevelance
num_columns <- ncol(abu_IBD_filtered) 

Prev_IBD <- numeric(nrow(abu_IBD_filtered))
nrow(abu_IBD_filtered)

for (i in 1:nrow(abu_IBD_filtered)){
  
  Prev_IBD[i] <- 100-((zero_counts[i]/num_columns)*100)
}



df_prev_IBD <- data.frame(Species=rownames(abu_IBD_filtered),Prevalence=Prev_IBD[1:7])
df_prev_IBD <- df_prev_IBD %>% mutate(category = "Healthy adults (vegan diet)")




####### abundance plot ########
abu_combined_all <- bind_rows(
  abu_western_y_long %>% mutate(category = "Healthy adults (vegetarian diet)"),
  abu_western_n_long %>% mutate(category = "Healthy adults (omnivore diet)"),
  abu_base_long %>% mutate(category = "Healthy adults (all)"),
  abu_IBD_long %>% mutate(category = "Healthy adults (vegan diet)"))




# Defining breaks
my_breaks <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100)

# Defining labels 
my_labels <- c( "0,0001","0.0001", "0.001", "0.01", "0.1", "1","10","100")

ggplot(abu_combined_all, aes(x = Species, y = Abundance, fill = category)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), color = "black", size = 1, outlier.size = 2, outlier.shape=16, outlier.colour="black") +
  theme_minimal() +
  labs(title = "Relative abundance of species in stool samples (different habitual diets)",
       x = "Species", y = "Relative abundance (%)") +
  theme(
    plot.title = element_text(hjust = 0.5, size=32),  
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 30),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28, color="black",face = "italic"),
    axis.text.y = element_text( size = 28, color="black"), 
    legend.position = "bottom", 
    legend.box = "horizontal", 
    legend.title = element_blank(),
    legend.text = element_text(size=27),
    legend.key.size = unit(2.5, "cm"),
    panel.grid.major.y = element_line(color = "gray65", size = 0.75),
    panel.grid.minor.y = element_line(color = "gray75", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  scale_y_log10(breaks = my_breaks, labels = my_labels) +
  scale_x_discrete(labels = species_names) 



####### Prevalence plot ##########
combined_all <- rbind(df_west_y,df_west_n,df_prev,df_prev_IBD)




ggplot(combined_all, aes(x = Species, y = Prevalence, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Prevalence of species in stool samples (different habitual diets)",
    x = "Species", y = "Prevalence (%)", fill = "Non-westernized"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size=34),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 32, color="black",face = "italic"),
    axis.text.y = element_text(size = 32, color="black"),
    legend.position = "bottom",
    legend.text = element_text(size=30),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),
    panel.grid.major.y = element_line(color = "gray65", size = 0.75),
    panel.grid.minor.y = element_line(color = "gray75", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  scale_x_discrete(labels = species_names) +
  scale_y_sqrt(
    breaks = seq(0, 100, by = 10),  
    limits = c(0, 100),            
    expand = c(0, 0) )








################ signifacns test #################


combined_data <- bind_rows(
  abu_western_y_long %>% mutate(category = "Healthy adults (vegetarian diet)"),
  abu_western_n_long %>% mutate(category = "Healthy adults (omnivore diet)"),
  abu_IBD_long %>% mutate(category = "Healthy adults (vegan diet)"))


#  Kruskal-Wallis test
kruskal_test_result <- combined_data %>%
  group_by(Species) %>%
  kruskal_test(Abundance ~ category)

# Filter Kruskal-Wallis test results by p-value
significant_species <- kruskal_test_result %>%
  filter(p < 0.05) %>%
  pull(Species)

# Subset data to include only significant phyla
filtered_data <- combined_data %>%
  filter(Species %in% significant_species)

# Dunn test only on significant phyla
dunn_test_result2 <- filtered_data %>%
  group_by(Species) %>%
  dunn_test(Abundance ~ category, p.adjust.method = "bonferroni")


write.table(dunn_test_result2, file = 'dunn_species.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
write.table(kruskal_test_result, file = 'kruskal_species.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")



