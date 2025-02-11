library(ggplot2)
library(dplyr)

# Define the T-number classification
hmmer_jr <- hmmer_jr %>%
  mutate(
    T_number = case_when(
      contig_length >= 2000 & contig_length < 5000  ~ "T=1",
      contig_length >= 5000 & contig_length < 8000  ~ "T≈1, T≈2",
      contig_length >= 8000 & contig_length < 12000 ~ "T=2",
      contig_length >= 12000 & contig_length < 14000 ~ "T≈3",
      contig_length >= 14000 & contig_length < 16000 ~ "T=3",
      contig_length >= 16000 & contig_length < 20000 ~ "T=4",
      TRUE ~ "Unknown"
    )
  )

# Check distribution of assigned T-numbers
table(hmmer_jr$T_number)

# Define T-number color scheme
t_colors <- c(
  "T=1" = "black",
  "T≈1, T≈2" = "purple",
  "T=2" = "red",
  "T≈3" = "orange",
  "T=3" = "gold",
  "T=4" = "yellow",
  "Unknown" = "gray"
)

# Create the histogram
ggplot(hmmer_jr, aes(x = contig_length, fill = T_number)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.8) +
  scale_fill_manual(values = t_colors) +
  labs(
    title = "Histogram of Contig Lengths with Assigned T-Numbers",
    x = "Contig Length (bp)",
    y = "Count",
    fill = "T-Number"
  ) +
  theme_minimal()
