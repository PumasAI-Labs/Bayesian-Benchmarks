library(arrow)
library(readr)

convert_csv_to_arrow <- function(csv_files, arrow_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(arrow_dir)) {
    dir.create(arrow_dir)
  }
  
  # Loop through CSV files
  for (csv_file in csv_files) {
    # Read CSV file
    df <- read_csv(csv_file, comment = "#")
    
    # Define Arrow file name and path
    arrow_file <- file.path(arrow_dir, paste0(basename(csv_file), ".arrow"))
    
    # Write Arrow file
    arrow::write_feather(df, arrow_file)
  }
}

