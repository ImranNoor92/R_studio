########################################################################################
hmmer <- read.csv("hmmer_clean_by_t.csv")
#Hmmer data from hmmer_clean_by_t.csv is upto date the cleanest, neatest data.
# I will be using these dataset for further analysis (02/18/2025)




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
library(stringr)
library(dplyr)

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

# View first few rows of the updated dataset
head(hmmer_function)

# Save as CSV
write.csv(hmmer, "jr_fold_classified.csv", row.names = FALSE)
jr_fold_classified <- read.csv("jr_fold_classified.csv")
