# Load required packages
library(dplyr)
library(readr)
# 1. Read and combine your search-result TSV files
# Adjust the pattern if your file names differ.
tsv_files <- list.files(pattern = "^SearchResults-.*\\.tsv$")

# Read each file and combine them into one data frame.
search_results <- lapply(tsv_files, function(file) {
  read_tsv(file, col_types = cols())
}) %>% 
  bind_rows()

# 2. Remove duplicate entries based on the Accession column.
search_results <- search_results %>% 
  distinct(Accession, .keep_all = TRUE)

# 3. Filter to include only entries from PFAM in the Source Database column.
search_results <- search_results %>% 
  filter(`Source Database` == "PFAM")
# Rename the 'Accession' column to lowercase 'accession'
# so that it matches the column name in your reference file.
search_results <- search_results %>% 
  rename(accession = Accession)
search_results <- search_results %>%
  rename(name = Name, pfam_description = Description)
# 3. Read the additional PFAM entries from the CSV file.
old_pfam <- read.csv("jr_pfam.csv", header = TRUE, stringsAsFactors = FALSE)
# Rename the column containing PFAM accessions to "accession" for consistency.
old_pfam <- old_pfam %>% 
  rename(accession = PFAM.Accession.Number)
old_pfam <- old_pfam %>%
  rename(name = Domain.Name)
write.csv(search_results, "search_results.csv", row.names = FALSE)
write.csv(old_pfam, "old_pfam.csv", row.names = FALSE)
# 4. Combine PFAM accession values from both sources.
#    We use union() to get the unique PFAM accession values from both datasets.
all_pfams_df <- full_join(search_results, old_pfam, by = "accession")
all_pfams <- union(search_results$accession, old_pfam$accession)
all_pfams_df <- data.frame(accession = all_pfams, stringsAsFactors = FALSE)
 HEAD
all_pfams_df <- read.csv("pfam_list_detailed.csv")
colnames(all_pfams_df)

levels(all_pfams_df$accession)
all_pfams_df$accession <- as.factor(all_pfams_df$accession)

# Or, if you prefer using tibble (from the tidyverse):
library(tibble)
all_pfams_df <- tibble(accession = all_pfams)
# 5. Read your main reference file.
#    Replace "your_reference_file.tsv" with the actual filename.
reference <- read.csv("All_DTRs_domtblout.csv")
colnames(reference)
# 6. Clean the reference file’s 'accession' column by removing the decimal part.
#    For example, this turns "PF00001.1" into "PF00001".
reference <- reference %>%
  mutate(accession = sub("\\..*", "", accession))

# 7. Filter the reference file so that only rows with a PFAM accession in 'all_pfams' are retained.
final_reference <- reference %>% 
  filter(accession %in% all_pfams_df$accession)
#make acession column in final_reference factor then check levels
final_reference$accession <- as.factor(final_reference$accession)
final_reference$target_name <- as.factor(final_reference$target_name)
levels(final_reference$accession)
write.csv(data.frame(accession = levels(final_reference$accession)), "unique_pfam.csv", row.names = FALSE)

# 8. Export the final filtered reference file as a CSV.
write_csv(final_reference, "final_joined.csv")
library(dplyr)

# Suppose your final_reference data frame has a column "query_name"
# Create a new column "query_name_clean" by removing an underscore and trailing digits:
final_reference <- final_reference %>%
  mutate(query_name_clean = sub("_[0-9]+$", "", query_name))

# Check the transformation:
head(final_reference$query_name_clean)
biome_info <- read.csv("DTRs_20kb.csv")
# Now, perform the left join using the cleaned column from final_reference 
# and matching it to the "genome_id" column in biome_info.
hmmer <- left_join(final_reference, biome_info, by = c("query_name_clean" = "genome_id"))
write_csv(hmmer, "hmmer.csv")
# Check the structure of the resulting data frame:
str(hmmer)
