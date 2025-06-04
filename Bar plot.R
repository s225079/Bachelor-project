# Load necessary libraries



library(ggplot2)
library(reshape2)



data <- read.delim("C:/Users/Miras/OneDrive - Danmarks Tekniske Universitet/Skrivebord/6. semester/Bachelor/Genome analysis/Barplot.txt")


# Melt the data frame
data_melted <- melt(data, id.vars = 'Strain')
data_melted$variable <- gsub("\\.", "-", data_melted$variable)



substrate_colors <- c(
  'alpha-glucan' = 'lightpink',
  'arabinan'='orchid2',
  'beta-galactan'='olivedrab3',
  'beta-glucan' = 'mediumpurple1',
  'beta-mannan'='royalblue',
  'fructan'='orange',
  'pectin' = 'lightcoral',
  'starch' = 'turquoise',
  'xylan' = 'seagreen3'
)





#  stacked bar chart
ggplot(data_melted, aes(x = Strain, y = value, fill = variable)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 9) +
  scale_fill_manual(values = substrate_colors) + 
  labs(title = 'Counts of PULs for each substrate in bacterial strains', x = 'Bacterial strain', y = 'Total number of PULs', fill='Substrate') +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 28),
    plot.title = element_text(size = 30, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic",color="black"),
    axis.text.y= element_text(color="black"),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26),
    panel.grid.major.y = element_line(color = "grey80", size = 0.8)
    
  )


