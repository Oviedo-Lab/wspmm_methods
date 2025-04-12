
# Load library to source Python code snippets
if ( !is.element("package:reticulate", search()) ) { # only try to load if not already loaded
  if(!require(reticulate)) {                         # if not installed, install
    install.packages("reticulate")                   # install if needed
    library(reticulate)                              # load
  }
}

# What version of python are we using? (For grabbing ROI masks)
# ... use_python("/usr/local/bin/python3.12")
use_virtualenv("~/envs/allen-env", required = TRUE)

# Orientation of reference-atlas space: 
# ... * Anterior -> Posterior * Superior -> Inferior * Left -> Right

# Set paths to files
data_path <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/Allen_data.csv"
data_path_h5 <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/C57BL6J-638850-raw.h5ad"
data_path_cellmeta <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/cell_metadata.csv"

make_slice_data <- function(
    slice = 21, 
    gene.list = c("Bcl11b", "Rorb") 
  ) {
    
    # Load sparse-matrix data
    file <- hdf5r::H5File$new(data_path_h5, mode = "r")
    #list.groups(file)
    #list.datasets(file)
    
    # Grab data from H5 file
    gene_codes <- file[["/X"]][["indices"]][]
    cell_id_idx <- file[["/X"]][["indptr"]][]
    cell_id_list <- file[["/obs/cell_label"]][]
    slice_id_list <- file[["/obs/brain_section_label/categories"]][]
    slice_num_list <- file[["/obs/brain_section_label/codes"]][]
    trscrpt_id_list <- file[["/var/transcript_identifier"]][]
    trscrpt_gene_symb_list <- file[["/var/gene_symbol"]][]
    trscrpt_gene_id_list <- file[["/var/gene_identifier"]][]
    run_lengths <- diff(cell_id_idx)
    gene_idx <- gene_codes + 1
    
    # Mask gene_idx for this slice
    gene_idx_mask = rep(slice_num_list, times = run_lengths) == slice
    gene_idx <- gene_idx[gene_idx_mask]
    
    # Mask cell_id for this slice
    slice_mask <- slice_num_list == slice
    slice_num_list <- slice_num_list[slice_mask]
    cell_id_list <- cell_id_list[slice_mask]
    run_lengths <- run_lengths[slice_mask]
    
    # Reorganize masked data into data frame
    num_initial_rows <- sum(gene_idx_mask)
    data <- data.frame(
      trscrpt_ct = file[["/X"]][["data"]][gene_idx_mask],
      trscrpt_gene_symb = trscrpt_gene_symb_list[gene_idx],
      slice_num = rep(slice_num_list, times = run_lengths),
      region = rep(NA, num_initial_rows),
      cell_id = rep(cell_id_list, times = run_lengths),
      slice_id = rep(slice_id_list[slice_num_list + 1], times = run_lengths),
      trscrpt_id = trscrpt_id_list[gene_idx],
      trscrpt_gene_id = trscrpt_gene_id_list[gene_idx]
    )
    
    # Remove blanks
    blank_mask <- grepl("Blank", data$trscrpt_gene_symb)
    data <- data[!blank_mask, ]
    
    # Get cell locations and remove unlabeled cells
    cell_meta <- data.table::fread(data_path_cellmeta)
    locatable_cells <- unique(cell_meta$cell_label)
    data <- data[data$cell_id %in% locatable_cells, ]
    
    # Keep only genes of interest
    data <- data[data$trscrpt_gene_symb %in% gene.list, ]
    
    # Put spatial coordinates into data frame
    data <- merge(
      data,
      cell_meta[, c("cell_label", "x", "y", "z")],
      by.x = "cell_id",
      by.y = "cell_label",
      all.x = TRUE
    )
    
    # Reorganize columns
    data <- data[,c("trscrpt_ct", "trscrpt_gene_symb", "slice_num", "region", "x", "y", "z", "cell_id", "slice_id", "trscrpt_id", "trscrpt_gene_id")]
    
    # Return data 
    return(data)
    
  }

slice_data_back <- make_slice_data(slice = 30)
slice_data_mid <- make_slice_data(slice = 35)
slice_data_front <- make_slice_data(slice = 40)

# genes_as_int <- as.integer(factor(slice_data$trscrpt_gene_symb))
# gene_colors <- rep("gray", nrow(slice_data)) 
# gene_colors[genes_as_int == 1] <- "red"
# gene_colors[genes_as_int == 2] <- "blue"
#plot(jitter(slice_data$x), jitter(slice_data$y), col = gene_colors, pch = 19, cex = 0.5)

# These dimensions go x, z, y? 
source_python("Allen_CCF.py")
masks <- generate_roi_masks()
dim_names <- c("AntPost", "SupInf", "LefRig")
for (m in 1:length(masks)) {
  colnames(masks[[m]]) <- dim_names
}

# slice_ROI_mask <- masks$ROI_mask_S1_L23[,2] == 416
# slice_ROI_mask <- TRUE
# ROI <- masks$ROI_mask_S1_L23[slice_ROI_mask, c(1,3)]
# points(ROI[,2]/100, ROI[,1]/100, col = "black", pch = 19, cex = 0.5)

library(rgl)
mult <- 100
plot3d(
  x = masks$ROI_mask_S1_L23[,1],
  y = masks$ROI_mask_S1_L23[,2] + 200,
  z = masks$ROI_mask_S1_L23[,3],
  col = "gray",
  ylim = c(0, 1000), xlim = c(0, 1000), zlim = c(0, 1000),
)
plot3d(
  x = (13.4 - slice_data_back$z) * mult, 
  y = slice_data_back$y * mult,
  z = slice_data_back$x * mult,
  col = "red",
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (13.4 - slice_data_mid$z) * mult, 
  y = slice_data_mid$y * mult,
  z = slice_data_mid$x * mult,
  col = "blue",
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (13.4 - slice_data_front$z) * mult, 
  y = slice_data_front$y * mult,
  z = slice_data_front$x * mult,
  col = "green",
  add = TRUE,
  aspect = TRUE
)


