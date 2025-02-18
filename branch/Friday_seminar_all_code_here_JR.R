########################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(viridis)
library(stringr)
library(dplyr)
library(patchwork)
#hmmer <- read.csv("hmmer_clean_by_t.csv")
#Hmmer data from hmmer_clean_by_t.csv is upto date the cleanest, neatest data.
# I will be using these dataset for further analysis (02/18/2025)

hmmer <- read.csv("jr_fold_classified.csv")
#Cleanest data so far


########################################################################################

# Method 1: Using genes_viral > 0 (Strict Gene-Based Filtering)
# Filter dataset where genes_viral > 0
hmmer$origin_genes <- ifelse(hmmer$genes_viral > 0, "viral", "non_viral")
# Do not use this method. Too much inclusive for viral genes.
# Create a dataset based on gene presence
jr_viral_genes <- subset(hmmer, origin_genes == "viral")
jr_non_viral_genes <- subset(hmmer, origin_genes == "non_viral")

# Summarize counts
viral_summary_genes <- table(hmmer$origin_genes)

# Convert to dataframe for visualization
viral_summary_genes_df <- as.data.frame(viral_summary_genes)

# View summary
print(viral_summary_genes_df)
#Var1 Freq
#1 non_viral    6
#2     viral  706
# View first few rows
head(jr_viral_genes)
head(jr_non_viral_genes)
#
#
########################################################################################
# Define keywords for viral origin in description column
# Step 1: 
# Load required libraries

# Convert all descriptions to lowercase, trim spaces, and fix encoding
hmmer$description_clean <- tolower(trimws(iconv(hmmer$description, to="ASCII//TRANSLIT")))

# Fix potential formatting issues (remove extra spaces, handle hyphens)
hmmer$description_clean <- str_replace_all(hmmer$description_clean, "-", " ")
hmmer$description_clean <- str_replace_all(hmmer$description_clean, "\\(", "\\\\(")  # Escape (
hmmer$description_clean <- str_replace_all(hmmer$description_clean, "\\)", "\\\\)")  # Escape )

# Define refined keyword lists with improved formatting
viral_keywords <- c("\\bvirus\\b", "\\bviral\\b", "\\bbacteriophage\\b", "\\bcapsid\\b", 
                    "\\bcoat protein\\b", "\\bprophage\\b", "\\bpicornavirus\\b", 
                    "\\bnepovirus\\b", "\\bwaikavirus\\b", "\\bpolyhedrin\\b", 
                    "\\bhexon\\b", "\\bvp4\\b")

prokaryotic_keywords <- c("\\bbcla c terminal domain\\b", "\\bpeptide n glycosidase f, c terminal\\b", 
                          "\\bendospore appendages\\b", "\\bcsp protease b prodomain\\b", 
                          "\\bclostridium neurotoxin, n terminal receptor binding\\b", 
                          "\\bexotoxin a binding\\b", "\\byersinia pseudotuberculosis mitogen\\b", 
                          "\\bbacteria\\b", "\\bbacterial\\b", "\\bprokaryotic\\b", 
                          "\\bribosomal\\b", "\\boperon\\b", "\\bpeptide n glycosidase\\b", 
                          "\\bcsp protease\\b", "\\bcarbohydrate esterase\\b", 
                          "\\bhomogentisate\\b", "\\byersinia\\b", "\\bclostridium\\b", "\\bendospore\\b")

eukaryotic_keywords <- c("\\bzpr1 zinc finger domain\\b", "\\bcip1 like, core domain\\b", 
                         "\\baspartyl asparaginyl beta hydroxylase\\b", "\\breelin subrepeat b\\b", 
                         "\\btgf beta propeptide\\b", "\\bpira like\\b")

prok_euk_viral_origin <- c("\\bcarbohydrate binding module 64\\b", "\\bthioredoxin like c terminal domain\\b", 
                           "\\bcarbohydrate binding module \\(family 35\\)\\b")  # Parentheses fixed

non_viral_origin <- c("\\beutq like cupin domain\\b", "\\bpolysaccharide lyase\\b", 
                      "\\bhomogentisate 1,2 dioxygenase n terminal\\b", "\\blyase, n terminal\\b", 
                      "\\balpha l arabinofuranosidase b, catalytic\\b", "\\bcarbohydrate esterase 2 n terminal\\b", 
                      "\\bcollagen binding domain\\b", "\\bgalactose binding domain like\\b", 
                      "\\bl proline 3 hydroxylase, c terminal\\b", 
                      "\\bphosphomannose isomerase type i, catalytic domain\\b", 
                      "\\bglycosyl hydrolase family 65, c terminal domain\\b", "\\bpeptidase family a21\\b")

cellular_keywords <- c("\\bribosome\\b", "\\bmembrane bound\\b", "\\bstructural component\\b")

# Ensure all keyword lists are lowercase
viral_keywords <- tolower(viral_keywords)
prokaryotic_keywords <- tolower(prokaryotic_keywords)
eukaryotic_keywords <- tolower(eukaryotic_keywords)
prok_euk_viral_origin <- tolower(prok_euk_viral_origin)
non_viral_origin <- tolower(non_viral_origin)
cellular_keywords <- tolower(cellular_keywords)

# Assign origin categories with correct prioritization
hmmer$origin <- case_when(
  str_detect(hmmer$description_clean, paste(viral_keywords, collapse = "|")) ~ "Viral",
  str_detect(hmmer$description_clean, paste(prokaryotic_keywords, collapse = "|")) ~ "Prokaryotic",
  str_detect(hmmer$description_clean, paste(eukaryotic_keywords, collapse = "|")) ~ "Eukaryotic",
  str_detect(hmmer$description_clean, paste(cellular_keywords, collapse = "|")) ~ "Cellular",
  str_detect(hmmer$description_clean, paste(prok_euk_viral_origin, collapse = "|")) ~ "Varies (Prok/Euk/Viral)",
  str_detect(hmmer$description_clean, paste(non_viral_origin, collapse = "|")) ~ "Varies (Prok/Euk)",
  TRUE ~ "Unknown"
)

# Debugging: Check why certain keywords are not being detected
hmmer %>%
  filter(str_detect(description_clean, paste(prok_euk_viral_origin, collapse = "|"))) %>%
  select(description_clean, origin) %>%
  head(100)

# Create dataset with assigned origins
hmmer_origin <- hmmer

# Summarize counts of origins
origin_summary <- table(hmmer_origin$origin)

# Convert summary to a dataframe
origin_summary_df <- as.data.frame(origin_summary)

# Print summary
print(origin_summary_df)

# View first few rows
head(hmmer_origin)
hmmer_origin <- write.csv(hmmer_origin,"hmmer_origin.csv")

########################################################################################

# Convert descriptions to lowercase, trim spaces, and normalize encoding
hmmer$description_clean <- tolower(trimws(iconv(hmmer$description, to="ASCII//TRANSLIT")))

# Standardize description formatting (fix hyphens, parentheses)
hmmer$description_clean <- str_replace_all(hmmer$description_clean, "-", " ")  # Convert hyphen to space
hmmer$description_clean <- str_replace_all(hmmer$description_clean, "\\(", "\\\\(")  # Escape (
hmmer$description_clean <- str_replace_all(hmmer$description_clean, "\\)", "\\\\)")  # Escape )

# Define refined keyword lists
capsid_keywords <- c("\\bcapsid\\b", "\\bcoat protein\\b", "\\bvirion\\b", "\\bbacteriophage tail\\b", 
                     "\\bviral coat\\b", "\\bhexon\\b", "\\bpolyhedrin\\b", "\\bvp4\\b", 
                     "\\bmajor structural protein\\b")

enzyme_keywords <- c("\\bhydrolase\\b", "\\bpolymerase\\b", "\\bkinase\\b", "\\boxidoreductase\\b", 
                     "\\blyase\\b", "\\bpeptidase\\b", "\\bphosphomannose isomerase\\b", 
                     "\\bglycosyl hydrolase\\b", "\\baspartyl asparaginyl\\b", 
                     "\\bpeptide n glycosidase\\b", "\\bcarbohydrate esterase\\b", 
                     "\\bhomogentisate\\b", "\\base\\b", "\\bpolysaccharide lyase\\b",
                     "\\bhomogentisate 1,2 dioxygenase n terminal\\b", "\\blyase, n terminal\\b", 
                     "\\balpha l arabinofuranosidase b, catalytic\\b", 
                     "\\bcarbohydrate esterase 2 n terminal\\b")

binding_keywords <- c("\\bbinding\\b", "\\bligand\\b", "\\breceptor\\b", "\\bcarbohydrate binding\\b", 
                      "\\bcollagen binding\\b", "\\bzinc finger\\b", "\\btoxin binding\\b", 
                      "\\bgalactose binding\\b", "\\bprotein binding\\b")

structural_keywords <- c("\\bmembrane\\b", "\\bcytoskeleton\\b", "\\bscaffolding\\b", 
                         "\\bnucleoplasmin\\b", "\\bendospore appendages\\b", 
                         "\\bstructural protein\\b", "\\bjelly roll\\b")

signaling_keywords <- c("\\bsignaling\\b", "\\bkinase\\b", "\\bphosphorylation\\b", 
                        "\\bresponse regulator\\b", "\\btgf beta\\b", "\\breelin\\b")

toxin_keywords <- c("\\btoxin\\b", "\\bexotoxin\\b", "\\bclostridium neurotoxin\\b", 
                    "\\bvirulence factor\\b", "\\bimmune modulation\\b")

unknown_keywords <- c("\\bduf\\b", "\\bunknown function\\b", "\\bpira like\\b")

# Convert lists to lowercase
capsid_keywords <- tolower(capsid_keywords)
enzyme_keywords <- tolower(enzyme_keywords)
binding_keywords <- tolower(binding_keywords)
structural_keywords <- tolower(structural_keywords)
signaling_keywords <- tolower(signaling_keywords)
toxin_keywords <- tolower(toxin_keywords)
unknown_keywords <- tolower(unknown_keywords)

# Assign function categories with prioritization
hmmer$functional_category <- case_when(
  str_detect(hmmer$description_clean, regex(paste(capsid_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Capsid Protein",
  str_detect(hmmer$description_clean, regex(paste(enzyme_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Enzyme",
  str_detect(hmmer$description_clean, regex(paste(binding_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Binding Protein",
  str_detect(hmmer$description_clean, regex(paste(structural_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Structural Protein",
  str_detect(hmmer$description_clean, regex(paste(signaling_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Signaling Protein",
  str_detect(hmmer$description_clean, regex(paste(toxin_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Toxin/Virulence",
  str_detect(hmmer$description_clean, regex(paste(unknown_keywords, collapse = "|"), ignore_case = TRUE)) ~ "Unknown",
  TRUE ~ "Other"
)

# Debugging: Check if function keywords are being detected correctly
hmmer %>%
  filter(str_detect(description_clean, paste(capsid_keywords, collapse = "|"))) %>%
  select(description_clean, functional_category) %>%
  head(20)

# Create a dataset with assigned functions
hmmer_function <- hmmer

# Summarize function categories
function_summary <- table(hmmer_function$functional_category)

# Convert summary to a dataframe
function_summary_df <- as.data.frame(function_summary)

# Print summary
print(function_summary_df)
write.csv(hmmer_function,"hmmer_function.csv")
# View first few rows of the updated dataset
head(hmmer_function)
#hmmer_function <- write.csv(hmmer_function,"hmmer_function.csv")


########################################################################################
#Plotting
#load all data
jrf_classified <- read.csv("jr_fold_classified.csv")
#Plot the origin
# Count occurrences of each origin category
origin_summary_df <- jrf_classified %>%
  count(origin)  # Creates a count dataframe

# Create the correct bar plot
origin_plot <- ggplot(origin_summary_df, aes(x = n, y = fct_reorder(origin, n), fill = origin)) +
  geom_col() +  # Use geom_col() to plot precomputed counts
  scale_fill_viridis_d(option = "magma") +
  theme_minimal() +
  labs(title = "Distribution of JRF Origins",
       x = "Count", y = "Origin") +
  theme(legend.position = "none")
origin_plot
#Plot the function
# Count occurrences of each functional category
function_summary_df <- jrf_classified %>%
  count(functional_category)  # Creates a count dataframe

# Create the correct bar plot
function_plot <-ggplot(function_summary_df, aes(x = n, y = fct_reorder(functional_category, n), fill = functional_category)) +
  geom_col() +  # Use geom_col() to plot precomputed counts
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal() +
  labs(title = "Distribution of Functional Categories in JRFs",
       x = "Count", y = "Functional Category") +
  theme(legend.position = "none")
function_plot
#Plot the origin and function together
# Arrange the plots side by side
combined_plot <- origin_plot + function_plot + plot_layout(ncol = 2)

# Display the combined plot
combined_plot

########################################################################################
#Viral plotting
jr_viral <- read.csv("jrf_viral.csv")
colnames(jr_viral)
# Summarize counts for each Biome1, Biome2, and T-number
# Summarize counts for each Biome1, Biome2, and T-number
heatmap_data <- jr_viral %>%
  count(biome1, biome2, T_number)  # Keep T_number as a variable
heatmap_data$T_number <- as.factor(heatmap_data$T_number)

heatmap_data$short_labels <- substr(interaction(heatmap_data$biome1, heatmap_data$biome2, sep = " - "), 1, 30)

ggplot(heatmap_data, aes(x = T_number, y = short_labels)) +
  geom_tile(aes(fill = n), color = "white") +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal(base_size = 14) +
  labs(title = "Heatmap of Biome1, Biome2, and T-number",
       x = "T-number",
       y = "Biome Categories",
       fill = "Count") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )
  ##########################################################################################
# Create the heatmap with Biome1 on X-axis, Biome2 on Y-axis, and T-number as color
ggplot(heatmap_data, aes(x = biome1, y = biome2, fill = T_number)) +
  geom_tile(color = "white") +  # Add white grid lines for clarity
  scale_fill_viridis_d(option = "plasma") +  # Plasma color for softer contrast
  theme_minimal(base_size = 14) +
  labs(title = "Biome1 vs Biome2 with T-number as Color",
       x = "Biome 1",
       y = "Biome 2",
       fill = "T-number") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate X-axis labels
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

heatmap_data$T_number <- as.numeric(as.character(heatmap_data$T_number))

# Create the heatmap with a smooth transition
ggplot(heatmap_data, aes(x = biome1, y = biome2, fill = T_number)) +
  geom_tile(color = NA) +  # Remove grid lines for a smoother appearance
  scale_fill_gradientn(colors = c(low = "darkgray", high = "orange"), 
                       name = "T-number") +  # Soft purple gradient
  theme_minimal(base_size = 14) +
  labs(title = "Smooth Heatmap of Biome1 vs Biome2 with T-number",
       x = "Biome 1",
       y = "Biome 2",
       fill = "T-number") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "grey90", color = NA),  # Light grey background for soft blending
    panel.grid.major = element_blank(),  # Remove grid lines for a smoother effect
    panel.grid.minor = element_blank()
  )
  ########################################################################################


library(dplyr)

# Reassign T-number based on contig_length
jr_viral <- jr_viral %>%
  mutate(T_number = case_when(
    contig_length <= 11500 ~ 1,   # T=1 for contigs ≤ 11,500
    contig_length <= 15000 ~ 3,   # T=3 for contigs >11,500 and ≤ 15,000
    contig_length <= 20000 ~ 4,   # T=4 for contigs >15,000 and ≤ 20,000
    TRUE ~ NA_real_              # Assign NA if outside range
  ))

# Convert T_number to a factor for better visualization
jr_viral$T_number <- factor(jr_viral$T_number, levels = c(1, 3, 4))

# Check if reassignment is correct
table(jr_viral$T_number)

# View first few rows after reassignment
head(jr_viral)



jr_viral$accession <- as.factor(jr_viral$accession)
levels(jr_viral$accession)
str(jr_viral$accession)
multiple_t_box <- ggplot(jr_viral, aes(x = accession, y = as.factor(T_number))) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.4) +  # Transparent boxplot
  geom_jitter(width = 0.2, size = 3, aes(color = as.factor(T_number))) +  # Scatter effect
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title = "Variation of T-numbers within Pfam Accessions",
       x = "Pfam Accession",
       y = "T-number",
       color = "T-number") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
multiple_t_box
########################################################################################

# Generate synthetic data
set.seed(42)
contig_lengths <- sample(500:20000, 500, replace=TRUE)

# Create DataFrame
df <- data.frame("Contig_Length_bp" = contig_lengths)

# Assign T-number
df <- df %>%
  mutate(T_number = case_when(
    Contig_Length_bp <= 11500 ~ "T=1",
    Contig_Length_bp <= 15000 ~ "T=3",
    Contig_Length_bp <= 20000 ~ "T=4",
    TRUE ~ "Unknown"
  ))

# Plot histogram with bins of 500
ggplot(df, aes(x = Contig_Length_bp, fill = T_number)) +
  geom_histogram(binwidth = 500, alpha=0.7, color="black") +
  scale_fill_viridis_d(option = "magma") +
  labs(title = "Distribution of Contig Lengths with T-number Assignments",
       x = "Contig Length (bp)", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################################################################################
# Load required libraries
library(dplyr)
library(ggplot2)

# Identify Pfam accessions that have both T=1 and either T=3 or T=4
pfam_multiple_t <- jr_viral %>%
  group_by(accession) %>%
  summarize(unique_t = list(unique(T_number))) %>%
  filter(any(c(1,3) %in% unique_t) | any(c(1,4) %in% unique_t)) %>%
  pull(accession)

# Filter the dataset to only include those Pfam accessions
jr_viral_filtered <- jr_viral %>%
  filter(accession %in% pfam_multiple_t)

# View the first few rows to confirm
head(jr_viral_filtered)
# Plot Biome1 vs Biome2, colored by T-number
ggplot(jr_viral_filtered, aes(x = biome1, fill = as.factor(T_number))) +
  geom_bar(position = "dodge") +
  scale_fill_viridis_d(name = "T-number") +
  theme_minimal() +
  labs(title = "Distribution of Pfam Accessions with Multiple T-numbers by Biome",
       x = "Biome1", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ########################################################################################
ggplot(jr_viral_filtered, aes(x = biome1, fill = as.factor(T_number))) +
  geom_bar(position = "dodge") +
  scale_fill_viridis_d(name = "T-number") +
  theme_minimal() +
  labs(title = "T-number Distribution Across Biomes for Multi-T Pfam Accessions",
       x = "Biome1", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~accession, scales = "free_y")  # Separate plots for each Pfam Accession
#####################################
library(ggplot2)
library(dplyr)

# Filter dataset to only include Pfam accessions that appear in multiple T-numbers
pfam_t_counts <- jr_viral_filtered %>%
  group_by(accession) %>%
  summarise(unique_t = n_distinct(T_number)) %>%
  filter(unique_t > 1)

# Subset original dataset to include only these accessions
multi_t_pfam <- jr_viral_filtered %>%
  filter(accession %in% pfam_t_counts$accession)

# Scatter plot of T-number, faceted by Pfam accession
library(ggplot2)
library(dplyr)

# Filter dataset to only include Pfam accessions that appear in multiple T-numbers
pfam_t_counts <- jr_viral_filtered %>%
  group_by(accession) %>%
  summarise(unique_t = n_distinct(T_number)) %>%
  filter(unique_t > 1)

# Subset original dataset to include only these accessions
multi_t_pfam <- jr_viral_filtered %>%
  filter(accession %in% pfam_t_counts$accession)

# Scatter plot of T-number, faceted by Pfam accession
ggplot(multi_t_pfam, aes(x = as.factor(T_number), y = biome1, color = as.factor(T_number))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d(name = "T-number") +
  theme_minimal() +
  labs(title = "Pfam Accessions with Multiple T-numbers",
       x = "T-number", y = "Biome1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~accession, scales = "free_y")  # Facet by Pfam accession
