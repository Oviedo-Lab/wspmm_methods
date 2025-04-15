
# Load library to source Python code snippets
if ( !is.element("package:reticulate", search()) ) { # only try to load if not already loaded
  if(!require(reticulate)) {                         # if not installed, install
    install.packages("reticulate")                   # install if needed
    library(reticulate)                              # load
  }
}

# Load rgl for 3D plotting
library(rgl)

# What version of python are we using? (For grabbing ROI masks)
# ... use_python("/usr/local/bin/python3.12")
use_virtualenv("~/envs/allen-env", required = TRUE)

# Orientation of reference-atlas space: 
# ... * Anterior -> Posterior * Superior -> Inferior * Left -> Right

# Set paths to files
data_path <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/Allen_data.csv"
data_path_h5 <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/C57BL6J-638850-raw.h5ad"
data_path_cellmeta <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/cell_metadata.csv"

# Helper function to get slice data
make_slice_data <- function(
    slice = 21, 
    gene.list = c("Bcl11b", "Rorb", "Gad2", "Foxb1", "Pvalb", "Slc38a1") 
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
      hemi = rep(NA, num_initial_rows),
      layer = rep(NA, num_initial_rows),
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
    # ... random downsampling, for when doing registration
    # n_rows <- nrow(data)
    # test_rows <- sample(c(1:n_rows), size = as.integer(n_rows/10), replace = FALSE)
    # data <- data[test_rows, ]
    
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
    data <- data[,c("trscrpt_ct", "trscrpt_gene_symb", "slice_num", "hemi", "layer", "x", "y", "z", "cluster_id", "cell_id", "slice_id", "trscrpt_id", "trscrpt_gene_id")]
    
    # Return data 
    return(data)
    
  }

# Make slice data
slice_data <- list()
slice_range <- c(30, 31, 33, 35, 36, 37, 38)
# ... 30:40 cover S1, but throw out 32, 34, 39, 40 for bad quality
for (s in seq_along(slice_range)) {
  cat("\nSlice data, slice", slice_range[s])
  slice_data[[s]] <- make_slice_data(slice = slice_range[s])
}
names(slice_data) <- paste0("slice_", slice_range)

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

# Define helper functions
logistic1 <- function(x, a = 10) {
  return((2/(1 + exp(-a*x))) - 1)
}
arctan <- function(x, y) {
  # x = rise, y = run
  theta <- atan(abs(y/x))
  theta[x < 0 & y > 0] <- pi - theta[x < 0 & y > 0]
  theta[x < 0 & y < 0] <- pi + theta[x < 0 & y < 0]
  theta[x > 0 & y < 0] <- 2*pi - theta[x > 0 & y < 0]
  return(theta)
}
line_slopeInter <- function(origin, point) {
  intercept <- (point["x"]*origin["y"] - point["y"]*origin["x"])/(point["x"] - origin["x"])
  slope <- (point["y"] - intercept)/point["x"]
  param <- c(slope, intercept)
  names(param) <- c("slope", "intercept")
  return(param)
}
find_midpoint <- function(p1, p2) {
  names(p1) <- c("x", "y") 
  names(p2) <- c("x", "y")
  x_ <- p1["x"] + (p2["x"] - p1["x"])/2 
  y_ <- p1["y"] + (p2["y"] - p1["y"])/2
  p <- c(x_, y_)
  names(p) <- c("x", "y")
  return(p)
}
find_topbot <- function(this_side, this_layer, base_mask, superinf_adj) {
  SupInf_ <- c(masks$SupInf[base_mask])[this_side & this_layer] + superinf_adj
  LefRig_ <- c(masks$LefRig[base_mask])[this_side & this_layer]
  Top_idx <- which.max(SupInf_)
  Top_ <- c(
    LefRig_[Top_idx], 
    SupInf_[Top_idx]
  )
  Bot_idx <- which.min(SupInf_)
  Bot_ <- c(
    LefRig_[Bot_idx], 
    SupInf_[Bot_idx]
  )
  out <- cbind(Top_, Bot_)
  colnames(out) <- c("Tp", "Bt")
  rownames(out) <- c("x", "y")
  out <- as.data.frame(out)
  return(out)
}
vec_dist <- function(point, point_set) {
  point <- array(point, dim = c(1,length(point)))
  colnames(point_set) <- c("x", "y")
  colnames(point) <- c("x", "y")
  dist <- sqrt((unlist(point[, "x"]) - point_set[, "x"])^2 + (unlist(point[, "y"]) - point_set[,"y"])^2)
  return(dist)
}

# Define geometric transformations
to_polar <- function(p) {
  if (is.null(dim(p))) {
    p <- array(p, dim = c(1,length(p)))
  }
  colnames(p) <- c("x", "y")
  r <- sqrt(p$x^2 + p$y^2)
  theta <- arctan(p$x, p$y)
  return(data.frame(r, theta))
}
to_cart <- function(p) {
  if (is.null(dim(p))) {
    p <- array(p, dim = c(1,length(p)))
  }
  colnames(p) <- c("r", "theta")
  x <- p$r * cos(p$theta)
  y <- p$r * sin(p$theta)
  return(data.frame(x, y))
}
translate <- function(p, p0, p1 = NULL, center = c(0, 0)) {
  colnames(p) <- c("x", "y")
  names(p0) <- c("x", "y")
  if (is.null(p1)) {
    p$x <- p$x - p0["x"]
    p$y <- p$y - p0["y"]
  } else {
    names(p1) <- c("x", "y")
    names(center) <- c("x", "y")
    mask_left <- p$x < center["x"]
    mask_right <- !mask_left
    p$x[mask_left] <- p$x[mask_left] - p0["x"]
    p$y[mask_left] <- p$y[mask_left] - p0["y"]
    p$x[mask_right] <- p$x[mask_right] - p1["x"]
    p$y[mask_right] <- p$y[mask_right] - p1["y"]
  }
  return(p)
}
bend <- function(p, p0, p0L, p0R, sL, sR, AL, AR) {
  
  # Scale the bend 
  max_value <- max(p) 
  AL <- AL/max_value 
  AR <- AR/max_value

  # Translate to origin and L/R centers of mass
  p_ <- translate(p, p0)
  pL <- translate(p, p0L)
  pR <- translate(p, p0R)
  
  # Find angles to rotate coordinates about CoM to align y-axis with given slope
  thetaL <- arctan(c(sL), c(1))
  thetaR <- arctan(c(sR), c(1))
  
  # Rotate coordinates to align y-axis with given slope
  pL <- to_polar(pL)
  pR <- to_polar(pR)
  pL$theta <- pL$theta + thetaL
  pR$theta <- pR$theta + thetaR
  pL <- to_cart(pL)
  pR <- to_cart(pR)
  
  # Bend coordinates around y-axis
  pL$y <- pL$y + AL * pL$x^2
  pR$y <- pR$y - AR * pR$x^2
  
  # Undo rotation 
  pL <- to_polar(pL)
  pR <- to_polar(pR)
  pL$theta <- pL$theta - thetaL
  pR$theta <- pR$theta - thetaR
  pL <- to_cart(pL)
  pR <- to_cart(pR)
  
  # Translate each hemisphere back to the origin
  mask_left <- p_$x < 0
  p[mask_left,] <- translate(pL, -p0L)[mask_left,]
  p[!mask_left,] <- translate(pR, -p0R)[!mask_left,]
  
  return(p)
  
}
rotate <- function(p, p0, p0L, p0R, AL, AR) {
  
  # Translate to origin and L/R centers of mass
  p_ <- translate(p, p0)
  pL <- translate(p, p0L)
  pR <- translate(p, p0R)
  
  # Convert to polar coordinates
  pLp <- to_polar(pL)
  pRp <- to_polar(pR)
  
  # Perform rotation
  pLr <- data.frame(r = pLp$r, theta = pLp$theta + AL)
  pRr <- data.frame(r = pRp$r, theta = pRp$theta + AR)
  
  # Convert back to cartesian coordinates
  pL <- to_cart(pLr)
  pR <- to_cart(pRr)

  # Translate each hemisphere back to the origin
  mask_left <- p_$x < 0
  p[mask_left,] <- translate(pL, -p0L)[mask_left,]
  p[!mask_left,] <- translate(pR, -p0R)[!mask_left,]
  
  return(p)
}
pinch <- function(p, p0, gammaL, gammaR) {
  p <- translate(p, p0)
  pp <- to_polar(p)
  ppr2 <- pp$r^2
  pL <- data.frame(r = sqrt(ppr2 * gammaL), theta = pp$theta)
  pR <- data.frame(r = sqrt(ppr2 * gammaR), theta = pp$theta)
  pL <- to_cart(pL)
  pR <- to_cart(pR)
  max_right <- max(p$x)
  max_left <- min(p$x)
  mask_left <- p$x < 0
  mask_right <- !mask_left
  left_ratio <- p$x[mask_left] / max_left
  right_ratio <- p$x[mask_right] / max_right
  p$x[mask_left] <- pL$x[mask_left]
  p$y[mask_left] <- pL$y[mask_left]
  p$x[mask_right] <- pR$x[mask_right]
  p$y[mask_right] <- pR$y[mask_right]
  p <- translate(p, -p0)
  return(p)
}

# Define function to transform and plot slices
process_slice <- function(
    slicedata, 
    transL = c(0, 0),
    transR = c(0, 0),
    pinchL = 1,
    pinchR = 1,
    rotateL = 0,
    rotateR = 0,
    bendL = 0,
    bendR = 0,
    make_plots = FALSE,
    threeD = FALSE,
    slice_num = NULL
  ) {
    
    # ... For registration, use "1/20" and randomly sample transcript rows for representative image of slice
    downsample_ratio <- 1/10
    
    # Find the z coordinates of these slices 
    z_coord <- rep(0, length(slicedata))
    
    # Initialize 3D plot with masks
    if (!make_plots) threeD <- FALSE
    if (threeD) {
      these_rows <- sample(c(1:nrow(masks)), size = as.integer(nrow(masks)*downsample_ratio/2), replace = FALSE)
      plot3d(
        x = masks$AntPost[these_rows],
        y = c(masks$SupInf + 180)[these_rows],
        z = masks$LefRig[these_rows],
        col = "steelblue",
        ylim = c(0, 1000), xlim = c(0, 1000), zlim = c(0, 1000),
      )
    }
    
    # Process each coronal slice
    for (s in seq_along(slicedata)) {
      
      # ... skip if not indicated
      if (!is.null(slice_num)) {
        if (slice_range[s] != slice_num) { next}
      }
      
      if (make_plots) cat("\nPlotting slice", slice_range[s])
      else cat("\nProcessing slice", slice_range[s])
      
      # Down-sample data
      if (make_plots) {
        slicedata_s <- slicedata[[s]][
          seq(
            1, nrow(slicedata[[s]]), 
            length.out = as.integer(nrow(slicedata[[s]])*downsample_ratio)
          ), ]
      } else {
        slicedata_s <- slicedata[[s]]
      }
      
      # Convert slice coordinates into CCF
      mult <- 100
      xm <- 1320/mult
      ym <- 800/mult
      zm <- 1140/mult
      slicedata_s$x <- (zm - slicedata_s$x) * mult
      slicedata_s$y <- (slicedata_s$y) * mult
      slicedata_s$z <- (xm - slicedata_s$z) * mult
      
      # Grab the z coordinate of the slice
      z_coord[s] <- as.integer(slicedata_s$z[1])
      
      # Plot slice 
      if (threeD) {
        plot3d(
          x = slicedata_s$z, 
          y = slicedata_s$y,
          z = slicedata_s$x,
          col = "gray",
          add = TRUE,
          aspect = TRUE
        )
      } else {
        
        # Grab z and set down-sample size
        mask_z <- masks$AntPost == z_coord[s]
        ds_size <- 1000 
        superinf_adj <- 180
        
        # Find center of mass of whole slice and label hemispherse
        center <- c(
          min(slicedata_s$x) + (max(slicedata_s$x) - min(slicedata_s$x))/2, 
          min(slicedata_s$y) + (max(slicedata_s$y) - min(slicedata_s$y))/2
        )
        left_rows <- slicedata_s$x < center[1]
        slicedata_s[left_rows, "hemi"] <- "left"
        slicedata_s[!left_rows, "hemi"] <- "right"
        
        # Find the center of mass of the mask in each hemisphere
        Lmask_mask <- masks$AntPost == z_coord[s] & masks$LefRig < center[1]
        Rmask_mask <- masks$AntPost == z_coord[s] & masks$LefRig > center[1]
        center_Lmask <- c(
          min(masks$LefRig[Lmask_mask]) + (max(masks$LefRig[Lmask_mask]) - min(masks$LefRig[Lmask_mask]))/2, 
          min(masks$SupInf[Lmask_mask]) + (max(masks$SupInf[Lmask_mask]) - min(masks$SupInf[Lmask_mask]))/2
        )
        center_Rmask <- c(
          min(masks$LefRig[Rmask_mask]) + (max(masks$LefRig[Rmask_mask]) - min(masks$LefRig[Rmask_mask]))/2, 
          min(masks$SupInf[Rmask_mask]) + (max(masks$SupInf[Rmask_mask]) - min(masks$SupInf[Rmask_mask]))/2
        )
        
        # Find slopes of Left and Right regions of interest
        left_side <- masks[mask_z, "LefRig"] < center[1]
        layer1 <- masks[mask_z, "Layer"] == "L1"
        layer6 <- masks[mask_z, "Layer"] == "L6a"
        L1TB <- find_topbot(left_side, layer1, mask_z, superinf_adj)
        L6TB <- find_topbot(left_side, layer6, mask_z, superinf_adj)
        R1TB <- find_topbot(!left_side, layer1, mask_z, superinf_adj)
        R6TB <- find_topbot(!left_side, layer6, mask_z, superinf_adj)
        L1mid <- find_midpoint(L1TB$Tp, L1TB$Bt)
        L6mid <- find_midpoint(L6TB$Tp, L6TB$Bt)
        R1mid <- find_midpoint(R1TB$Tp, R1TB$Bt)
        R6mid <- find_midpoint(R6TB$Tp, R6TB$Bt)
        Lnorm_slopeInter <- line_slopeInter(L1mid, L6mid)
        Rnorm_slopeInter <- line_slopeInter(R1mid, R6mid)
        Lslope <- Lnorm_slopeInter["slope"]
        Rslope <- Rnorm_slopeInter["slope"]
        
        # Coordinate Transforms for alignment with CCF ("registration")
        slicedata_s[,c("x", "y")] <- translate(slicedata_s[,c("x", "y")], transL, transR, center)
        slicedata_s[,c("x", "y")] <- rotate(slicedata_s[,c("x", "y")], center, center_Lmask, center_Rmask, rotateL, rotateR)
        slicedata_s[,c("x", "y")] <- bend(slicedata_s[,c("x", "y")], center, center_Lmask, center_Rmask, Lslope, Rslope, bendL, bendR)
        slicedata_s[,c("x", "y")] <- pinch(slicedata_s[,c("x", "y")], center, pinchL, pinchR)
        
        # Plot background slice
        if (make_plots) {
          plot_title <- paste0("Slice ", slice_range[s], ", z = ", z_coord[s])
          plot(slicedata_s$x, slicedata_s$y, col = "gray", pch = 19, cex = 0.5, asp = 1, main = plot_title, ylim = c(0, 1000), xlim = c(0, 1200))
        }
        
        # Find candidate points for annotation
        leftright <- masks[mask_z, "LefRig"]
        superinf <- masks[mask_z, "SupInf"] + superinf_adj
        min_x <- min(leftright)
        max_x <- max(leftright)
        min_y <- min(superinf)
        max_y <- max(superinf)
        candidate_rows <- slicedata_s$x < max_x & slicedata_s$x > min_x & slicedata_s$y < max_y & slicedata_s$y > min_y
        slicedata_s <- slicedata_s[candidate_rows, ]
        layer_distances <- array(NA, dim = c(nrow(slicedata_s), length(Layer_names)))
        
        # Extract layers from mask
        if (make_plots) ds_array <- array(0, dim = c(ds_size*length(Layer_names), 3))
        ntrans <- nrow(slicedata_s)
        cat("\nROI search,", ntrans, "transcripts: ")
        tracker <- as.integer(seq(1, ntrans, length.out = 10))
        for (i in seq_along(Layer_names)) {
          
          # Grab coordinates in this layer and z slice
          slidemask <- masks$Layer == Layer_names[i] & mask_z
          if (sum(slidemask) == 0) next
          
          # Subset mask to slice
          leftright <- masks[slidemask, "LefRig"]
          superinf <- masks[slidemask, "SupInf"] + superinf_adj
          
          
          # threshold = 0.05
          
          # Identify transcripts in data falling near these points 
          cat("\nLayer", Layer_names[i], " ")
          mask_df <- data.frame(x = leftright, y = superinf)
          tracker <- as.integer(seq(1, ntrans, length.out = 10))
          for (j in seq_along(1:ntrans)) {
            if (any(j == tracker)) cat("*")
            these_dist <- vec_dist(slicedata_s[j, c("x", "y")], mask_df)
            these_dist <- these_dist[these_dist != 0]
            layer_distances[j, i] <- -min(these_dist)
          }
          
          # Down-sample 
          n_slidemask <- sum(slidemask)
          these_rows <- sample(c(1:n_slidemask), size = ds_size, replace = n_slidemask < ds_size)
          row_range <- ((i-1)*ds_size + 1):(i*ds_size)
          
          # Fill array
          if (make_plots) {
            ds_array[row_range, 1] <- leftright[these_rows]
            ds_array[row_range, 2] <- superinf[these_rows]
            ds_array[row_range, 3] <- i + 1
          }
          
        }
        
        # Extract layer annotations 
        row_mins <- apply(-layer_distances, 1, FUN = min)
        inside_ROI <- row_mins < 0.05 * mult
        slicedata_s$layer <- Layer_names[max.col(layer_distances)]
        
        # Make plots or save data
        if (make_plots) {
          points(ds_array[,1], ds_array[,2], col = ds_array[,3], pch = 19, cex = 0.5)
          # Plot center point 
          names(center) <- c("x", "y")
          points(center["x"], center["y"], col = "black", pch = 19, cex = 3)
        } else {
          slicedata[[s]] <- slicedata_s[inside_ROI, ]
        }
        
      }
      
    }
    
    if (!make_plots) {
      return(slicedata)
    }
    
  }

# ... throw out, 32, 34, 39, 40

# Slice 30
slice_data <- process_slice(
  slice_data, 
  pinchL = 1.09, pinchR = 1.15, 
  rotateL = -0.075, rotateR = 0.1, 
  bendL = 0.05, bendR = 0.25,
  slice_num = 30
  )

# Slice 31
slice_data <- process_slice(
  slice_data, 
  transL = c(0, -20), transR = c(0, 0),
  pinchL = 1.01, pinchR = 0.96, 
  rotateL = -0.1, rotateR = 0.08, 
  bendL = -0.075, bendR = 0,
  slice_num = 31
)

# Slice 33
slice_data <- process_slice(
  slice_data, 
  transL = c(0, -45), transR = c(0, -10),
  pinchL = 1.025, pinchR = 0.95, 
  rotateL = -0.05, rotateR = 0.03, 
  slice_num = 33
)

# Slice 35
slice_data <- process_slice(
  slice_data, 
  transL = c(0, -68), transR = c(0, -32),
  pinchL = 0.97, pinchR = 1.02, 
  rotateL = -0.055, rotateR = 0.04, 
  bendL = -0.45, bendR = 0.1,
  slice_num = 35
)

# Slice 36
slice_data <- process_slice(
  slice_data, 
  pinchL = 0.95, pinchR = 1, 
  rotateL = -0.025, rotateR = 0.15, 
  bendL = 0, bendR = 0,
  slice_num = 36
)

# Slice 37
slice_data <- process_slice(
  slice_data, 
  transL = c(0, -60), transR = c(0, 0),
  pinchL = 1.01, pinchR = 0.95, 
  rotateL = -0.025, rotateR = 0.06, 
  bendL = 0, bendR = -0.35,
  slice_num = 37
)

# Slice 38 
slice_data <- process_slice(
  slice_data, 
  transL = c(0, -50), transR = c(0, -30),
  pinchL = 0.99, pinchR = 1, 
  rotateL = -0.06, rotateR = 0.16, 
  bendL = 0, bendR = 0,
  slice_num = 38
)




# our columnar axis is approx 1mm, these slices are 200um apart, should should be able to get 4 slices 
# Need to estimate the SupInf coordinate of our horizontal slices!!!??
























