
library(highcharter)
library(dplyr)
library(htmlwidgets)
library(webshot)
setwd("branch/jr")
dir()
hmmer_jr <- read.csv("hmmer_jr.csv")
# Aggregate counts for each Biome Level
summary_biome1 <- hmmer_jr %>% 
  group_by(biome1) %>% 
  summarise(size = n(), .groups = "drop")

summary_biome2 <- hmmer_jr %>% 
  group_by(biome1, biome2) %>% 
  summarise(size = n(), .groups = "drop")

summary_biome3 <- hmmer_jr %>% 
  group_by(biome1, biome2, biome3) %>% 
  summarise(size = n(), .groups = "drop")

# Convert Data to Hierarchical Format for Sunburst
dout_biome1 <- data_to_hierarchical(summary_biome1, c("biome1"), "size")
dout_biome2 <- data_to_hierarchical(summary_biome2, c("biome1", "biome2"), "size")
dout_biome3 <- data_to_hierarchical(summary_biome3, c("biome1", "biome2", "biome3"), "size")


pl_biome1 <- highchart() %>%
  hc_chart(type = "sunburst") %>%
  hc_title(text = "Biome Level 1 Sunburst Chart") %>%
  hc_add_series(name = "Biome Level 1", data = dout_biome1) %>%
  hc_plotOptions(series = list(dataLabels = list(format = "{point.name}"))) %>%
  hc_add_theme(hc_theme_economist())

pl_biome1


pl_biome2 <- highchart() %>%
  hc_chart(type = "sunburst") %>%
  hc_title(text = "Biome Level 2 Sunburst Chart") %>%
  hc_add_series(name = "Biome Level 2", data = dout_biome2) %>%
  hc_plotOptions(series = list(dataLabels = list(format = "{point.name}"))) %>%
  hc_add_theme(hc_theme_538())

pl_biome2

pl_biome3 <- highchart() %>%
  hc_chart(type = "sunburst") %>%
  hc_title(text = "Biome Level 3 Sunburst Chart") %>%
  hc_add_series(name = "Biome Level 3", data = dout_biome3) %>%
  hc_plotOptions(series = list(dataLabels = list(format = "{point.name}"))) %>%
  hc_add_theme(hc_theme_alone())
k
pl_biome3
