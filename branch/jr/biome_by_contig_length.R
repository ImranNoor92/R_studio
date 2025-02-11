library(ggplot2)
library(dplyr)

# Remove the "Unknown" category
biome_counts <- biome_counts %>%
  filter(biome1 != "Unknown") %>%  # Exclude Unknown category
  arrange(biome1, desc(count))  # Ensure sorting: Biome1 (A-Z) → Biome2 (by count)

# Convert Biome2 into a factor with correct sorting
biome_counts$biome2 <- factor(
  biome_counts$biome2, 
  levels = unique(biome_counts$biome2)
)

# Define color scheme (remove Unknown)
biome_colors <- c(
  "Aquatic" = "#1f78b4",
  "Terrestrial" = "#b15928",
  "Host_Associated" = "#33a02c",
  "Engineered" = "#e31a1c"
)

# Set x-axis range from 5000 to 20000 (Linear Scale)
x_min <- 0
x_max <- 20000

# Create the improved plot (No Log Scale, No Unknown)
ggplot(biome_counts, aes(y = biome2, x = avg_contig_length, fill = biome1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = biome_colors) +  
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set manual x-axis range
  geom_text(aes(label = count), hjust = -0.2, size = 5) +  
  theme_minimal() +
  labs(
    x = "Average Contig Length (bp)", 
    y = "Biome", 
    fill = "Biome Category",
    title = "Top Counts Across Main Biomes"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  )