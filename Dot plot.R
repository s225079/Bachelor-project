
library(readxl)
library(ggplot2)
library(tidyr)

library(dplyr)  # Gør sikker på dplyr er loaded

# Read the data from the file
data <- read.delim("C:/Users/Miras/OneDrive - Danmarks Tekniske Universitet/Skrivebord/6. semester/dot.txt")

# Print the column names to check them
print(colnames(data))

data_long <- data %>%
  pivot_longer(cols = -1, names_to = "Bacterial_Strain", values_to = "Presence") %>% filter(Presence == "+")  # Keep only rows where Presence is "+"




# plot
ggplot(data_long, aes(x = Bacterial_Strain, y = `SCFA_product`)) +
  geom_point(size = 6) + 
  theme_minimal() +
  labs(title = " SCFA Products by Bacterial Strains",
       x = "Bacterial Strain",
       y = "SCFA Product") +
  theme(axis.text = element_text(size = 26), 
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic",color="black"),
        axis.text.y = element_text(color="black"))
