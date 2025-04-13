
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
    data <- data[data$cell_id %in% cell_meta$cell_label, ]
    
    # Keep only genes of interest
    data <- data[data$trscrpt_gene_symb %in% gene.list, ]
    
    # Grab cluster ids 
    data$cluster_id <- cell_meta$cluster_alias[match(data$cell_id, cell_meta$cell_label)]
    
    # Put spatial coordinates into data frame
    data <- merge(
      data,
      cell_meta[, c("cell_label", "x", "y", "z")],
      by.x = "cell_id",
      by.y = "cell_label",
      all.x = TRUE
    )
    
    # Reorganize columns
    data <- data[,c("trscrpt_ct", "trscrpt_gene_symb", "slice_num", "region", "x", "y", "z", "cluster_id", "cell_id", "slice_id", "trscrpt_id", "trscrpt_gene_id")]
    
    # Return data 
    return(data)
    
  }

slice_data <- list()
slice_range <- 30:40
for (s in seq_along(slice_range)) {
  cat("\nSlice data, slice", s)
  slice_data[[s]] <- make_slice_data(slice = slice_range[s])
}

# These dimensions go x, z, y? 
make_masks <- FALSE 
Layer_names <- c("L1", "L23", "L4", "L5", "L6a", "L6b")
if (make_masks) {
  source_python("Allen_CCF.py")
  masks_list <- generate_roi_masks()
  masks_list <- masks
  dim_names <- c("AntPost", "SupInf", "LefRig", "Layer")
  masks <- data.frame()
  for (m in seq_along(masks_list)) {
    masks_list[[m]] <- as.data.frame(masks_list[[m]])
    masks_list[[m]]$Layer <- Layer_names[m]
    colnames(masks_list[[m]]) <- dim_names
    masks <- rbind(masks, masks_list[[m]])
  }
  rm(masks_list)
  write.csv(masks, file = "masks.csv", row.names = FALSE)
} else {
  masks <- read.csv("masks.csv")
}

mult <- 100
xm <- 1320/mult
ym <- 800/mult
zm <- 1140/mult

plot_slice <- function(slicedata) {
  
  # Find the z coordinates of these slices 
  z_coord <- rep(0, length(slicedata))
  
  for (s in seq_along(slicedata)) {
    
    # Convert slice coordinates into CCF
    slicedata[[s]]$x <- (zm - slicedata[[s]]$x) * mult
    slicedata[[s]]$y <- (slicedata[[s]]$y) * mult
    slicedata[[s]]$z <- (xm - slicedata[[s]]$z) * mult
    
    # Grab the z coordinate of the slice
    z_coord[s] <- as.integer(slicedata[[s]]$z[1])
    
    # Plot background slice 
    plot(slicedata[[s]]$x, slicedata[[s]]$y, col = "gray", pch = 19, cex = 0.5)
    
    mask_z <- masks$AntPost == z_coord[s]
    for (i in seq_along(Layer_names)) {
      # Grab coordinates in this layer and z slice
      slidemask <- masks$Layer == Layer_names[i] & mask_z
      # Plot
      points(masks[slidemask, LefRig], masks[slidemask, SupInf] + 180, col = i, pch = 19, cex = 0.5)
    }
    
    
  }
  
}







# our columnar axis is approx 1mm, these slices are 200um apart, should should be able to get 4 slices 
# Need to estimate the SupInf coordinate of our horizontal slices!!!??

library(rgl)
plot3d(
  x = masks$ROI_mask_S1_L1[,1],
  y = masks$ROI_mask_S1_L1[,2] + 180,
  z = masks$ROI_mask_S1_L1[,3],
  col = "gray",
  ylim = c(0, 1000), xlim = c(0, 1000), zlim = c(0, 1000),
)

plot3d(
  x = (xm - slice_data_back$z) * mult, 
  y = (slice_data_back$y) * mult,
  z = (zm - slice_data_back$x) * mult,
  col = "red",
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_mid1$z) * mult, 
  y = (slice_data_mid1$y) * mult,
  z = (zm - slice_data_mid1$x) * mult,
  col = "blue4",
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_mid2$z) * mult, 
  y = (slice_data_mid2$y) * mult,
  z = (zm - slice_data_mid2$x) * mult,
  col = "blue2",
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_mid3$z) * mult, 
  y = (slice_data_mid3$y) * mult,
  z = (zm - slice_data_mid3$x) * mult,
  col = "blue",
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_front$z) * mult, 
  y = (slice_data_front$y) * mult,
  z = (zm - slice_data_front$x) * mult,
  col = "green",
  add = TRUE,
  aspect = TRUE
)

sss <- 10

plot3d(
  x = (xm - slice_data_backS$z) * mult, 
  y = (slice_data_backS$y) * mult,
  z = (zm - slice_data_backS$x) * mult,
  col = "black", size = sss,
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_mid1S$z) * mult, 
  y = (slice_data_mid1S$y) * mult,
  z = (zm - slice_data_mid1S$x) * mult,
  col = "black", size = sss,
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_mid2S$z) * mult, 
  y = (slice_data_mid2S$y) * mult,
  z = (zm - slice_data_mid2S$x) * mult,
  col = "black", size = sss,
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_mid3S$z) * mult, 
  y = (slice_data_mid3S$y) * mult,
  z = (zm - slice_data_mid3S$x) * mult,
  col = "black", size = sss,
  add = TRUE,
  aspect = TRUE
)
plot3d(
  x = (xm - slice_data_frontS$z) * mult, 
  y = (slice_data_frontS$y) * mult,
  z = (zm - slice_data_frontS$x) * mult,
  col = "black", size = sss,
  add = TRUE,
  aspect = TRUE
)

