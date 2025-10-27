
# Setup ################################################################################################################

# Install and load SPARK
# download gfortran: https://cran.r-project.org/bin/macosx/tools/
if (!require("SPARK", quietly = TRUE)) {devtools::install_github('xzhoulab/SPARK')}
library(SPARK)

# SPARK-X ##############################################################################################################

# Grab variable levels
last_var_col <- which(colnames(count_data) == "y_bins")
mice <- unique(count_data$mouse)
hemi <- unique(count_data$hemisphere)
sample_names <- interaction(mice, hemi, sep = "_")

age <- c()
for (m in mice) {
  age <- c(age, unique(as.character(count_data[count_data$mouse == m, "age"])))
}
age <- as.factor(age)

# Analyze with SPARK-X #################################################################################################

make_SPARKX_data <- function(
    h, 
    m, 
    cdata,
    laminar = FALSE # collapse bin coordinates along x axis
) {
  
  # Which cells (rows) should be grabbed? 
  mask <- cdata$mouse == m & cdata$hemisphere == h
  
  # Which genes (columns) should be grabbed? 
  last_var_col <- which(colnames(cdata) == "y_bins")
  gene_idx <- (last_var_col+1):ncol(cdata)
  
  # Which columns hold the bin coordinates? 
  bin_idx <- (last_var_col - 1):last_var_col
  
  # Prune down cdata 
  cnames <- colnames(cdata)
  gene_names <- cnames[gene_idx]
  bin_names <- cnames[bin_idx]
  cdata <- cdata[mask, c(gene_idx, bin_idx)]
  colnames(cdata) <- c(gene_names, bin_names)
  cdata_genes_only <- cdata[, gene_names]
  
  # Extract "spot" coordinates
  dim2c <- cdata[, bin_names[2]]
  if (laminar) dim1c <- rep(1, length(dim2c))
  else dim1c <- cdata[, bin_names[1]]
  dim2 <- 1:max(dim2c)
  dim1 <- 1:max(dim1c)
  
  # Make the count matrix for SPARKX
  # sp_count <- matrix(0, nrow = length(gene_names), ncol = length(dim1) * length(dim2))
  # spot_names <- c()
  # rownames(sp_count) <- gene_names
  # for (y in dim2) {
  #   for (x in dim1) {
  #     mask <- dim1c == x & dim2c == y
  #     spot_names <- c(spot_names, paste0(x, "x", y))
  #     sp_count[, (y-1)*length(dim1) + x] <- colSums(cdata_genes_only[mask, ], na.rm = TRUE)
  #   }
  # }
  # colnames(sp_count) <- spot_names
  
  # Make the count matrix for SPARKX
  # ... varied that this code here matches the loops above
  # Generate all possible spot names (ensuring full coverage)
  grid_index <- expand.grid(dim1, dim2)
  spot_names <- paste0(grid_index[,1], "x", grid_index[,2])
  # Initialize count matrix with zeros
  sp_count <- matrix(0, nrow = length(gene_names), ncol = length(spot_names))
  rownames(sp_count) <- gene_names
  colnames(sp_count) <- spot_names
  # Identify actual spots present in cdata_genes_only
  existing_spots <- paste0(dim1c, "x", dim2c)
  # Compute column index for each existing spot
  spot_lookup <- match(existing_spots, spot_names)
  # Loop over genes to fill in counts correctly
  for (i in seq_along(gene_names)) {
    counts <- tapply(cdata_genes_only[, i], existing_spots, sum, na.rm = TRUE)
    sp_count[i, match(names(counts), spot_names)] <- counts
  }
  
  # Make the location matrix for SPARKX
  info <- cbind.data.frame(
    x=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",1)),
    y=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",2))
  )
  rownames(info) <- colnames(sp_count)
  location <- as.matrix(info)
  
  # Return the remade matrices 
  return(list(sp_count = sp_count, location = location))
  
}

run_SPARKX <- function(
    h, 
    m, 
    cdata,
    laminar = FALSE
) {
  
  sdata <- make_SPARKX_data(h, m, cdata, laminar)
  sp_count <- sdata$sp_count
  location <- sdata$location
  
  sparkX <- sparkx(
    sp_count,
    location,
    numCores = ncores,
    option = "mixture"
  )
  
  return(sparkX)
  
}

SPARKX_results <- list()
for (r in 1:nrow(coldata)) {
  h <- coldata[r, "hemi"]
  m <- coldata[r, "mouse"]
  SPARKX_results[[paste0(h, "_", m)]] <- run_SPARKX(h, m, count_data)
}
