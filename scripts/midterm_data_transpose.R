library(dplyr)
library(tidyr)
library(readr)

raw_dat <- readr::read_csv("raw_data/midterm_cleaned_spectra.csv")

# Get column names from your actual data
cols <- colnames(raw_dat)

# Get sample names (columns that don't start with "...")
sample_names <- cols[!grepl("^\\.\\.\\.", cols)]

# First, collect all unique wavelengths and create sample data frames
sample_dfs <- list()

for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  
  # Find the column index of this sample
  sample_idx <- which(cols == sample_name)
  
  # The x-values are in the column before it
  x_col <- sample_idx - 1
  
  # Get x and y values
  x_vals <- raw_dat[[x_col]]
  y_vals <- raw_dat[[sample_idx]]
  
  # Create a data frame for this sample
  sample_dfs[[sample_name]] <- data.frame(
    wavelength = x_vals,
    value = y_vals
  ) %>%
    filter(!is.na(wavelength))  # Remove NA wavelengths
  
  colnames(sample_dfs[[sample_name]]) <- c("wavelength", sample_name)
}

# Get all unique wavelengths across all samples
all_wavelengths <- sort(unique(unlist(lapply(sample_dfs, function(x) x$wavelength))))

# Create base data frame with all wavelengths
result <- data.frame(wavelength = all_wavelengths)

# Join each sample to the result
for (sample_name in sample_names) {
  result <- result %>%
    left_join(sample_dfs[[sample_name]], by = "wavelength")
}

# Check result
dim(result)
head(result[, 1:5])
tail(result[, 1:5])


result_transposed <- result %>%
  tibble::column_to_rownames("wavelength") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample")

write.table(result_transposed, "processed_data/midterm_raw_data_transposed.csv", quote = F, row.names = F, col.names = T, sep = ",")






