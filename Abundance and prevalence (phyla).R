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
##### Collection of data #######

data <- sampleMetadata


dat <-sampleMetadata |>
  filter(age >= 18) |>
  filter(body_site == "stool") |>
  filter(disease == "healthy") |>
  select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance", rownames = "short")


dat1 <- colData(dat) %>% data.frame()

tax <- rowData(dat) %>% data.frame()

abu <- assay(dat) %>% data.frame()



n_base <- (nrow(dat1)/nrow(data))*100


######## Merga tax and abu dataframes #####



# Convert rowData (taxonomy) and assay (abundance) to data frames
tax <- rowData(dat) %>% data.frame()  # Extract taxonomy
abu <- assay(dat) %>% data.frame()    # Extract abundance

# Move row names of abu to a column to match taxonomy
abu <- abu %>% rownames_to_column(var = "species")

merged_data <- merge(tax[, c("species", "phylum")], abu, by = "species")



######## calculate abundance on phylum level #########

phylum_abundance <- merged_data %>%
  group_by(phylum) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) # Sum values per phylum

########## Plot ##########


phylum_abundance_long <- phylum_abundance %>%
  pivot_longer(cols = -phylum, names_to = "sample", values_to = "abundance")
phylum_abundance_long <- na.omit(phylum_abundance_long)



ggplot(phylum_abundance_long, aes(x = phylum, y = abundance, fill=phylum)) +
  geom_boxplot(outliers = FALSE) +
  theme_minimal() +
  labs(title = "Abundance of phyla in stool samples (healthy adults)",
       x = "Phyla", y = "Abundance (log10-transformed)" )+
  theme(
    plot.title = element_text(hjust = 0.5, size=22),  # Center the title
    axis.title.x = element_text(size = 20),  # Increase x-axis title font size
    axis.title.y = element_text(size = 20),  # Increase y-axis title font size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color="black"),
    axis.text.y = element_text(size = 20, color="black"),# Increase x-axis text size
    legend.position = "none") +
  scale_y_log10(labels = function(x) parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x))))  



########## Prevalence ###########



######### Plots of prevelance ########

#Counting the number of columns=0
zero_counts <- numeric(nrow(phylum_abundance))  # Create a vector to store results

for (i in 1:nrow(phylum_abundance)){
  zero_counts[i] <- sum(phylum_abundance[i, ] == 0)
}
#Calculating the prevelance
num_columns <- ncol(phylum_abundance)  # Number of columns

Prev <- numeric(nrow(phylum_abundance))

for (i in 1:nrow(phylum_abundance)){
  
  Prev[i] <- 100-((zero_counts[i]/num_columns)*100)
}


df_prev <- data.frame(phylum=phylum_abundance$phylum,Prevalence=Prev)

df_prev<- na.omit(df_prev)


############### Omnivore ##############



##### Collection of data #######




dat_IBD <-sampleMetadata |>
  filter(age >= 18) |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  filter(!is.na(diet)) |>
  filter(diet=="omnivore") |>
  returnSamples("relative_abundance", rownames = "short")


dat1_IBD <- colData(dat_IBD) %>% data.frame()



######## Merga tax and abu dataframes #####



# Convert rowData (taxonomy) and assay (abundance) to data frames
tax_IBD <- rowData(dat_IBD) %>% data.frame()  # Extract taxonomy
abu_IBD <- assay(dat_IBD) %>% data.frame()    # Extract abundance

# Move row names of abu to a column to match taxonomy
abu_IBD <- abu_IBD %>%
  rownames_to_column(var = "species")

merged_data_IBD <- merge(tax_IBD[, c("species", "phylum")], abu_IBD, by = "species")



######## calculate abundance on phylum level #########

phylum_abundance_IBD <- merged_data_IBD %>%
  group_by(phylum) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) # Sum values per phylum



# Convert to long format 
phylum_abundance_IBD_long <- phylum_abundance_IBD%>%
  pivot_longer(cols = -phylum, names_to = "sample", values_to = "abundance")

phylum_abundance_IBD_long <- na.omit(phylum_abundance_IBD_long)

########## Prevalence ###########


#Counting the number of columns=0
zero_counts <- numeric(nrow(phylum_abundance_IBD))  # Create a vector to store results

for (i in 1:nrow(phylum_abundance_IBD)){
  zero_counts[i] <- sum(phylum_abundance_IBD[i, ] == 0)
}


#Calculating the prevelance
num_columns <- ncol(phylum_abundance_IBD)  # Number of columns

Prev_IBD <- numeric(nrow(phylum_abundance_IBD))

for (i in 1:nrow(phylum_abundance_IBD)){
  
  Prev_IBD[i] <- 100-((zero_counts[i]/num_columns)*100)
}


df_prev_IBD <- data.frame(phylum=phylum_abundance_IBD$phylum,Prevalence=Prev_IBD)

df_prev_IBD <- na.omit(df_prev_IBD)

############### vegan ##############



##### Collection of data #######




dat_vegan <-sampleMetadata |>
  filter(age >= 18) |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  filter(!is.na(diet)) |>
  filter(diet=="vegan") |>
  returnSamples("relative_abundance", rownames = "short")


dat1_vegan <- colData(dat_vegan) %>% data.frame()



######## Merga tax and abu dataframes #####



# Convert rowData (taxonomy) and assay (abundance) to data frames
tax_vegan <- rowData(dat_vegan) %>% data.frame()  # Extract taxonomy
abu_vegan <- assay(dat_vegan) %>% data.frame()    # Extract abundance

# Move row names of abu to a column to match taxonomy
abu_vegan <- abu_vegan %>%
  rownames_to_column(var = "species")

merged_data_vegan <- merge(tax_vegan[, c("species", "phylum")], abu_vegan, by = "species")



########  abundance on phylum level #########

phylum_abundance_vegan <- merged_data_vegan %>%
  group_by(phylum) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) # Sum values per phylum



# Convert to long format
phylum_abundance_vegan_long <- phylum_abundance_vegan %>%
  pivot_longer(cols = -phylum, names_to = "sample", values_to = "abundance")

phylum_abundance_vegan_long <- na.omit(phylum_abundance_vegan_long)

########## Prevalence ###########



#Counting the number of columns=0
zero_counts <- numeric(nrow(phylum_abundance_vegan))  # Create a vector to store results

for (i in 1:nrow(phylum_abundance_vegan)){
  zero_counts[i] <- sum(phylum_abundance_vegan[i, ] == 0)
}


#Calculating the prevelance
num_columns <- ncol(phylum_abundance_vegan)  # Number of columns

Prev_vegan <- numeric(nrow(phylum_abundance_vegan))

for (i in 1:nrow(phylum_abundance_vegan)){
  
  Prev_vegan[i] <- 100-((zero_counts[i]/num_columns)*100)
}


df_prev_vegan <- data.frame(phylum=phylum_abundance_vegan$phylum,Prevalence=Prev_vegan)

df_prev_vegan <- na.omit(df_prev_vegan)




############### vegetarian ##############



##### Collection of data #######




dat_vege <-sampleMetadata |>
  filter(age >= 18) |>
  filter(body_site == "stool") |>
  select(where(~ !all(is.na(.x)))) |>
  filter(!is.na(diet)) |>
  filter(diet=="vegetarian") |>
  returnSamples("relative_abundance", rownames = "short")


dat1_vege <- colData(dat_vege) %>% data.frame()



######## Merga tax and abu dataframes #####



# Convert rowData (taxonomy) and assay (abundance) to data frames
tax_vege <- rowData(dat_vege) %>% data.frame()  # Extract taxonomy
abu_vege <- assay(dat_vege) %>% data.frame()    # Extract abundance

# Move row names of abu to a column to match taxonomy
abu_vege <- abu_vege %>%
  rownames_to_column(var = "species")

merged_data_vege <- merge(tax_vege[, c("species", "phylum")], abu_vege, by = "species")



######## calculate abundance on phylum level #########

phylum_abundance_vege <- merged_data_vege %>%
  group_by(phylum) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) # Sum values per phylum



# Convert to long format 
phylum_abundance_vege_long <- phylum_abundance_vege %>%
  pivot_longer(cols = -phylum, names_to = "sample", values_to = "abundance")

phylum_abundance_vege_long <- na.omit(phylum_abundance_vege_long)

########## Prevalence ###########



#Counting the number of columns=0
zero_counts <- numeric(nrow(phylum_abundance_vege))  # Create a vector to store results

for (i in 1:nrow(phylum_abundance_vege)){
  zero_counts[i] <- sum(phylum_abundance_vege[i, ] == 0)
}


#Calculating the prevelance
num_columns <- ncol(phylum_abundance_vege)  # Number of columns

Prev_vege <- numeric(nrow(phylum_abundance_vege))

for (i in 1:nrow(phylum_abundance_vege)){
  
  Prev_vege[i] <- 100-((zero_counts[i]/num_columns)*100)
}


df_prev_vege <- data.frame(phylum=phylum_abundance_vege$phylum,Prevalence=Prev_vege)

df_prev_vege <- na.omit(df_prev_vege)













############# combined ##########

abu_combined_all <- bind_rows(
  phylum_abundance_long %>% mutate(category = "Healthy adults (all)"),
  phylum_abundance_IBD_long %>% mutate(category = "Healthy adults (omnivorous diet)"),
  phylum_abundance_vegan_long %>% mutate(category = "Healthy adults (vegan diet)"),
  phylum_abundance_vege_long %>% mutate(category = "Healthy adults (vegetarian diet)"))

###### filter data##########



abu_combined_all <- abu_combined_all %>%
  filter(abundance != 0) %>%                             # remove zeroes
  group_by(phylum, category) %>%                         # group by phylum and category
  filter(n() > 1) %>%                                    # keep only groups with more than one datapoint
  ungroup()



# Define your breaks (manually or use log_breaks to generate them)
my_breaks <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100)

# Define labels with varying accuracy
my_labels <- c( "0,0001","0.0001", "0.001", "0.01", "0.1", "1","10","100")

########## new""""""""
ggplot(abu_combined_all, aes(x = phylum, y = abundance, fill = category)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), color = "black", size = 1, outlier.size = 2, outlier.shape=16, outlier.colour="black") +
  theme_minimal() +
  labs(title = "Relative abundance of phyla in stool samples (different habitual diets)",
       x = "Phyla", y = "Relative abundance (%)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 32),  # Center the title
    axis.title.x = element_text(size = 28),  # Increase x-axis title font size
    axis.title.y = element_text(size = 28),  # Increase y-axis title font size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28, color = "black",face = "italic"),
    axis.text.y = element_text(size = 28, color = "black"),  # Increase y-axis text size
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal",  # Arrange legend items horizontally
    legend.title = element_blank(),
    legend.text = element_text(size = 28),
    legend.key.size = unit(2, "cm"),
    panel.grid.major.y = element_line(color = "gray65", size = 0.75),
    panel.grid.minor.y = element_line(color = "gray75", size = 0.5),
    panel.grid.major.x = element_line(color="gray65", size=0.5)) +
  scale_y_log10(breaks = my_breaks, labels = my_labels) 


    
   

############# Prevalence plot ############





df_prev <- df_prev %>% mutate(category = "Healthy adults (all)")
df_prev_IBD <- df_prev_IBD %>% mutate(category = "Healthy adults (omnivorous diet)")
df_prev_vegan <- df_prev_vegan %>% mutate(category = "Healthy adults (vegan diet)")
df_prev_vege <- df_prev_vege %>% mutate(category = "Healthy adults (vegetarian diet)")

combined_all <- rbind( df_prev, df_prev_IBD, df_prev_vegan, df_prev_vege)


ggplot(combined_all, aes(x = phylum, y = Prevalence, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Prevalence of phyla in stool samples (different habitual diets)",
       x = "Phyla", y = "Prevalence (%)") +
  theme(
    plot.title = element_text(hjust = 0.5, size=32),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 30, color="black",face = "italic"),
    axis.text.y = element_text(size = 30, color="black"),
    legend.position = "bottom",
    legend.text = element_text(size=30),
    legend.title = element_blank() ,
    legend.key.size = unit(1.5, "cm"),
    panel.grid.major.y = element_line(color = "gray65", size = 0.75),
    panel.grid.minor.y = element_line(color = "gray75", size = 0.5)) +
  scale_y_sqrt(
    breaks = seq(0, 100, by = 10),  
    limits = c(0, 100),             
    expand = c(0, 0))



################ SIgnificans ###########

# Combine all long-form abundance data with diet label
phylum_abundance_IBD_long$diet <- "omnivore"
phylum_abundance_vegan_long$diet <- "vegan"
phylum_abundance_vege_long$diet <- "vegetarian"

combined_data <- bind_rows(
  phylum_abundance_IBD_long,
  phylum_abundance_vegan_long,
  phylum_abundance_vege_long
)



# Perform the Kruskal-Wallis test
kruskal_test_result <- combined_data %>%
  group_by(phylum) %>%
  kruskal_test(abundance ~ diet)

# Filter Kruskal-Wallis test results by p-value (e.g., p < 0.05)
significant_phyla <- kruskal_test_result %>%
  filter(p < 0.05) %>%
  pull(phylum)

# Subset the original data to include only significant phyla
filtered_data <- combined_data %>%
  filter(phylum %in% significant_phyla)

# Perform the Dunn test only on significant phyla
dunn_test_result2 <- filtered_data %>%
  group_by(phylum) %>%
  dunn_test(abundance ~ diet, p.adjust.method = "bonferroni")


write_xlsx(dunn_test_result2, "dunn_phylum.xlsx")
write_xlsx(kruskal_test_result, "kruskal_phylum.xlsx")


write.table(dunn_test_result2, file = 'dunn_phylum.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
write.table(kruskal_test_result, file = 'kruskal_phylum.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
