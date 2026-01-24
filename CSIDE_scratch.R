
# CSIDE benchmarking functions #########################################################################################

# Load library 
# options(timeout = 600000000) ### set this to avoid timeout error
# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)

# Function to convert data to CSIDE format and run RCTD
convert_sim_data_to_CSIDE <- function(
    sim,
    ref_counts,
    max_cores = 2
) {
  
  # Need:
  #  DATA:
  #   count matrix (one per "puck"): rows as genes, columns as spots (i.e., spatial locations, or "pixels")
  #   coords matrix: rows as spots, columns as x and y coordinates
  #   nUMI vector: total counts per spot
  #  REFERENCE:
  #   ref_counts: matrix with rows as genes, columns as cells, pseudo-bulk library
  #   cell_types: named vector with cell type for each cell in ref_counts
  #   nUMI vector: total counts per cell in ref_counts
  
  # Grab data and split replicates by fixed effect
  sim_data <- sim$data
  sim_data$replicate <- paste(sim_data$replicate, sim_data$fixedeffect, sep = "_")
  replicates <- sort(unique(sim_data$replicate))
  grounp_ids <- as.integer(grepl("_trt", replicates)) + 1 # 1 for ref, 2 for trt
  
  # full grid of pixels
  rx <- range(sim_data$bin_x)
  ry <- range(sim_data$bin_y)
  pixels <- CJ(
    bin_x = rx[1]:rx[2],
    bin_y = ry[1]:ry[2]
  )
  pixels <- pixels[, pixel := paste(bin_x, bin_y, sep = "_")]
  pixels <- as.matrix(pixels)
  
  # Convert pixel matrix into coords matrix 
  coords <- pixels 
  rownames(coords) <- coords[, "pixel"]
  coords <- as.data.frame(coords[, c("bin_x", "bin_y")])
  colnames(coords) <- c("x", "y")
  coords$x <- as.numeric(coords$x)
  coords$y <- as.numeric(coords$y)
  
  # Create cell types for REFERENCE
  cell_types <- data.frame(
    cell = colnames(ref_counts),
    type = sample(c("celltype1", "celltype2"), ncol(ref_counts), replace = TRUE)
  )
  cell_types <- setNames(cell_types[[2]], cell_types[[1]])
  cell_types <- as.factor(cell_types) # convert to factor data type
  
  # Randomly create ten filler genes that are differentially expressed between cell types
  n_filler_genes <- 10
  n_original_genes <- length(unique(sim_data$gene))
  ct_mask <- cell_types == "celltype1"
  n_ct1_ref <- sum(ct_mask)
  n_ct2_ref <- sum(!ct_mask)
  this_slice_num <- sim_data$slice_num[1]
  filler_df <- data.frame()
  original_gene_names <- rownames(ref_counts)
  for (fg in c(1:n_filler_genes)) {
    
    # Select a gene from ref_counts
    fg_name <- sample(original_gene_names, 1)
    ct1_counts <- ref_counts[fg_name, ct_mask]
    ct1_counts <- ct1_counts * runif(1) * 2
    de_mult <- (2 * runif(1))^2
    if (de_mult > 0.5 & de_mult <= 1) de_mult <- 0.5
    if (de_mult > 1 & de_mult <= 1.5) de_mult <- 1.5
    ct2_counts <- ct1_counts * de_mult
    fg_name <- paste0("_filler_", fg, "_", fg_name)
    
    # Update ref_counts with filler gene
    new_row <- rep(0, ncol(ref_counts))
    new_row[ct_mask] <- rpois(n_ct1_ref, ct1_counts + 1)
    new_row[!ct_mask] <- rpois(n_ct2_ref, ct2_counts + 1)
    ref_counts <- rbind(ref_counts, new_row)
    rownames(ref_counts)[nrow(ref_counts)] <- fg_name
    
    for (r in replicates) {
      replicate_num <- as.integer(sub(".*?(\\d+).*", "\\1", r))
      rep_rate_scalar <- sim$replicate_rate_scalars[replicate_num] + 0.5
      r_mask <- sim_data$replicate == r
      n_spots <- sum(r_mask)/n_original_genes
      r_idx <- sample(which(r_mask), n_spots, replace = TRUE) 
      r_df <- sim_data[r_idx,]
      ct_mask_r_df <- r_df$bin_x >= max(r_df$bin_x)/2 # arbitrarily assign left half to celltype1, right half to celltype2
      n_ct1 <- sum(ct_mask_r_df)
      n_ct2 <- n_spots - n_ct1
      r_df$count[ct_mask_r_df] <- rpois(n_ct1, sample(ct1_counts * 2 * rep_rate_scalar, n_ct1, replace = TRUE) + 1)
      r_df$count[!ct_mask_r_df] <- rpois(n_ct2, sample(ct2_counts * 2 * rep_rate_scalar, n_ct2, replace = TRUE) + 1)
      r_df$gene <- fg_name
      filler_df <- rbind(filler_df, r_df)
    }
    
  }
  sim_data <- rbind(sim_data, filler_df)
  
  # Make reference library counts 
  nUMI_ref <- colSums(ref_counts) 
  names(nUMI_ref) <- colnames(ref_counts)
  
  # Make reference
  reference <- Reference(ref_counts, cell_types, nUMI_ref)
  
  # Make data for each replicate
  n_replicates <- length(replicates)
  pucks <- list()
  length(pucks) <- n_replicates
  names(pucks) <- replicates
  barcodes <- list()
  length(barcodes) <- n_replicates
  names(barcodes) <- replicates
  for (r in replicates) {
    
    # Prune down data to just this replicate and necessary columns
    mask_r <- sim_data$replicate == r
    data_thin <- sim_data[mask_r,c("gene", "count", "bin_x", "bin_y")]
    # Count matrix 
    dt <- as.data.table(data_thin)
    # aggregate counts
    dt[, pixel := paste(bin_x, bin_y, sep = "_")]
    dt <- dt[, .(count = sum(count)), by = .(gene, pixel)]
    # ensure all gene–pixel combinations exist
    dt_full <- merge(
      CJ(gene = unique(dt$gene), pixel = pixels[,"pixel"]),
      dt,
      by = c("gene", "pixel"),
      all.x = TRUE
    )
    dt_full[is.na(count), count := 0]
    # cast to matrix
    pixel_mat <- dcast(
      dt_full,
      gene ~ pixel,
      value.var = "count",
      fill = 0
    )
    pixel_mat <- as.matrix(pixel_mat[, -1, with = FALSE])
    rownames(pixel_mat) <- dt_full[, unique(gene)]
    
    # Get library counts 
    nUMI_spatial <- colSums(pixel_mat)
    
    # Make and save puck
    pucks[[r]] <- SpatialRNA(coords, pixel_mat, nUMI_spatial)
    barcodes[[r]] <- colnames(pucks[[r]]@counts)
    
  }
  
  barcodes <- barcodes[[1]]
  
  # Create RCTD object
  low_panel_cutoff <- 0 # ... set to zero to account for small panel
  myRCTD.reps <- create.RCTD.replicates(
    spatialRNA.replicates = pucks,
    reference = reference,
    replicate_names = replicates,
    group_ids = grounp_ids,
    max_cores = max_cores,
    gene_cutoff = low_panel_cutoff, 
    fc_cutoff = low_panel_cutoff,
    gene_cutoff_reg = low_panel_cutoff, 
    fc_cutoff_reg = low_panel_cutoff,
    # UMI_min = low_panel_cutoff,
    # UMI_min_sigma = low_panel_cutoff,
    CELL_MIN_INSTANCE = 1,
    CONFIDENCE_THRESHOLD = 1
  )
  myRCTD.reps <- run.RCTD.replicates(myRCTD.reps)
  
  return(
    list(
      RCTD = myRCTD.reps,
      barcodes = barcodes
    )
  )
  
}

# Function to run CSIDE
run_CSIDE <- function(
    sim,
    ref_counts,
    max_cores = 2
) {
  
  # Convert data
  RCTD <- tryCatch(
    convert_sim_data_to_CSIDE(
      sim = sim,
      ref_counts = ref_counts,
      max_cores = max_cores
    ),
    error = function(e) {
      cat("\nError occurred with CSIDE data conversion or RCTD. Retrying with different random draws...")
      convert_sim_data_to_CSIDE(
        sim = sim,
        ref_counts = ref_counts,
        max_cores = max_cores
      )
    }
  )
  
  # Run CSIDE
  myRCTD.reps <- run.CSIDE.replicates(
    RCTD.replicates = RCTD$RCTD,
    cell_types = c("celltype1", "celltype2"),
    cell_type_threshold = 5, 
    weight_threshold = 0,
    fdr = 0.05,
    population_de = TRUE,
    de_mode = "nonparam",
    barcodes = RCTD$barcodes,
    log_fc_thresh = 0
  )
  
  # Run population inference with groups
  myRCTD.pop <- CSIDE.population.inference(
    myRCTD.reps, 
    fdr = 0.05, 
    use.groups = TRUE,
    MIN.CONV.REPLICATES = 1,
    MIN.CONV.GROUPS = 1,
    CT.PROP = 0,
    log_fc_thresh = 0
  )
  
  # Return completed model
  return(
    list(
      myRCTD.reps = myRCTD.reps,
      myRCTD.pop = myRCTD.pop
    )
  )
  
}

test <- run_CSIDE(
  sim = sim_data,
  ref_counts = ref_counts,
  max_cores = bs_chunksize
)

# CSIDE CODE

choose_sigma_gene <- function(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = F, n.iter = 100, MIN_CHANGE = 0.001, MAX_ITER_SIGMA = 10, PRECISION.THRESHOLD = .01) {
  sigma_s_best <- sigma_init
  sigma_vals <- names(Q_mat_all)
  alpha1 <- NULL; alpha2 <- NULL;
  MIN_ITERATIONS <- 15
  n.iter.tot <- 0
  for(iter in 1:MAX_ITER_SIGMA) {
    last_sigma <- sigma_s_best
    sigma_curr <- as.character(sigma_s_best)
    set_likelihood_vars(Q_mat_all[[sigma_curr]], X_vals, sigma = sigma_curr)
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter,
                                  MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                                  alpha1_init = alpha1, alpha2_init = alpha2, MIN_ITERATIONS = MIN_ITERATIONS)
    n.iter.tot <- n.iter.tot + res$n.iter
    alpha1 <- res$alpha1; alpha2 <- res$alpha2
    MIN_ITERATIONS <- 3
    prediction <- res$prediction
    pred_c <- as.vector(prediction)
    res_val <- numeric(length(sigma_vals))
    names(res_val) <- sigma_vals
    for(sigma_s in sigma_vals) {
      set_likelihood_vars(Q_mat_all[[as.character(sigma_s)]], X_vals, sigma = sigma_s)
      res_val[as.character(sigma_s)] <- (calc_log_l_vec_fast(pred_c, as.vector(t(Y))))
    }
    sigma_s_best <- names(which.min(res_val))
    if(sigma_s_best == last_sigma) {
      break
    }
  }
  res$n.iter <- n.iter.tot
  return(list(sigma_s_best = sigma_s_best, res = res))
}

mysweept <- function(tX2,tlk, K) {
  g_2 <- tX2[rep(1:dim(tX2)[1],K),] * tlk[rep(1:K, each = dim(tX2)[1]), ]
  return(g_2)
}

sweep1t <- function(tX1, lambda) {
  g_1 <- Rfast::eachrow(tX1, lambda,oper = "*")
  return(g_1)
}

sweep2t <- function(tX1, tdl, k) {
  if(dim(tX1)[1] > 0) {
    X1_Q <- Rfast::eachrow(tX1, tdl[k,],oper = "*")
  } else {
    X1_Q <- tX1
  }
  return(X1_Q)
}

sweep3t_all <- function(tX2, tdl, K) {
  X2_Q <- tX2[rep(1:dim(tX2)[1],K),] * tdl[rep(1:K, each = dim(tX2)[1]), ]
  return(X2_Q)
}

construct_hess_fast <- function(X1,X2,lambda,lambda_k, K, d1_d2) {
  tX1 <- t(X1); tX2 <- t(X2)
  g_1 <- sweep1t(tX1, lambda)
  tlk <- t(lambda_k)
  tdl <- Rfast::eachrow(tlk, d1_d2$d1_vec, '*')
  g_2 <- mysweept(tX2,tlk, K)
  grad <- rbind(g_1, g_2)
  grad_Q <- Rfast::eachrow(grad, d1_d2$d2_vec, '*')
  H1 <- -grad_Q %*% t(grad)
  X1_Q <- Rfast::eachrow(tX1, lambda*d1_d2$d1_vec, '*')
  grad_1 <- matrix(rowSums(X1_Q), 1, dim(X1)[2])
  H2_11 <- X1_Q %*% X1
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  H2_12 <- matrix(0,nrow = L1, ncol = L2*K)
  for(k in 1:K) {
    X1_Q <- sweep2t(tX1, tdl, k)
    H2_12[,(1 + L2*(k-1)):(L2*k)] <- X1_Q %*% X2
  }
  H2 <- matrix(0, nrow = L1 + L2*K, ncol = L1 + L2*K)
  if(L1 > 0) {
    H2[(1:L1),1:L1] <- H2_11
    H2[1:L1,(L1+1):(L1+L2*K)] <- H2_12
    H2[(L1+1):(L1+L2*K),1:L1] <- t(H2_12)
  }
  X2_Q <- sweep3t_all(tX2, tdl, K)
  grad_2 <- matrix(rowSums(X2_Q),dim(X2)[2],K)
  H2m <- X2_Q %*% X2
  for(k in 1:K) {
    H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <-
      H2m[(1 + (k-1)*L2):(k * L2),]# X2_Q %*% X2
  }
  H <- (H1 - H2)
  return(list(H = H, grad_1 = grad_1, grad_2 = grad_2))
}

# KEY FUNCTION 2
solveIRWLS.effects_trust <-function(Y, X1, X2, my_beta, test_mode, verbose = FALSE,
                                    n.iter = 200, MIN_CHANGE = .01, PRECISION.THRESHOLD = .05,
                                    alpha1_init = NULL, alpha2_init = NULL, MIN_ITERATIONS = 15){
  lam_threshold = 1e-8;
  beta_succ <- 1.1; beta_fail <- 0.5;
  gamma <- 0.8 # prev gamma = 0.1 / 0.25
  epsilon_2 <- 1e-6 * length(Y);
  delta <- 0.1 #trust region radius
  init_val <- min(-5, log(10/median(rowSums(my_beta))))
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  n_cell_types = dim(my_beta)[2]
  if(is.null(alpha1_init))
    alpha1 <- numeric(dim(X1)[2])
  else
    alpha1 <- alpha1_init
  if(is.null(alpha2_init)) {
    alpha2 <- matrix(0,nrow = dim(X2)[2], ncol = n_cell_types) # initialize it to be the previous cell type means
    alpha2[1,] <- init_val
    if(test_mode == 'categorical') {
      alpha2[,] <- init_val # multi mode
    }
  } else {
    alpha2 <- alpha2_init
  }
  pred_decrease_vals <- numeric(n.iter)
  Y[Y > K_val] <- K_val
  K <- dim(my_beta)[2]
  lambda_k <- exp(sweep(X2 %*% (alpha2), 1, X1 %*% (alpha1),'+'))*my_beta #J by K
  lambda <- rowSums(lambda_k)
  lambda[lambda < lam_threshold] <- lam_threshold
  d1_d2 <- calc_Q_all(Y, lambda)
  prev_ll <-  -sum(d1_d2$d0_vec)
  #prev_ll <- calc_log_l_vec(lambda,Y)
  error_vec <- (1:dim(my_beta)[2]) == 0
  for(itera in 1:n.iter) {
    H_list <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2)
    H <- H_list$H
    grad_1 <- H_list$grad_1; grad_2 <- H_list$grad_2
    d_vec_o <- c(grad_1, grad_2)
    D_mat_o <- psd(H)
    norm_factor <- norm(D_mat_o,"2")
    D_mat <- D_mat_o / norm_factor
    d_vec <- (d_vec_o / norm_factor) #- 2 * lambda_reg * c(alpha1, alpha2)
    #epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
    A <- cbind(diag(dim(D_mat)[2]), -diag(dim(D_mat)[2]))
    bzero <- rep(-delta,2*dim(D_mat)[2]); lambda_reg <- 1e-7
    D_mat <- D_mat + diag(dim(D_mat)[1])*lambda_reg # avoid numerical errors
    solution <-  quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution #* 0.01 # CHANGE THIS
    
    predicted_decrease = -(0.5*t(solution) %*% D_mat_o %*% solution - sum(d_vec_o*solution))
    
    alpha1_new <- alpha1 + solution[1:L1]
    alpha2_new <- alpha2 + matrix(solution[(L1 + 1):length(solution)], nrow = L2, ncol = K)
    lambda_k_new <- exp(sweep(X2 %*% (alpha2_new), 1, X1 %*% (alpha1_new),'+'))*my_beta #J by K
    error_vec <- is.na(colMeans(lambda_k_new)) | error_vec
    lambda_k_new[is.na(lambda_k_new)] <- 1
    lambda_new <- rowSums(lambda_k_new)
    d_all <- calc_Q_all(Y, lambda_new)
    #log_l <- calc_log_l_vec(lambda_new,Y)
    log_l <-  -sum(d_all$d0_vec)
    true_decrease = prev_ll - log_l
    
    pred_decrease_vals[itera] <- predicted_decrease
    if(true_decrease >= (gamma) * predicted_decrease) {
      delta = beta_succ*delta
      alpha1 <- alpha1_new
      alpha2 <- alpha2_new
      lambda_k <- lambda_k_new
      lambda <- lambda_new
      lambda[lambda < lam_threshold] <- lam_threshold
      prev_ll <- log_l
      d1_d2 <- d_all
      #d1_d2 <- get_d1_d2(Y, lambda)
    } else {
      delta = min(1,beta_fail*delta)
    }
    if(F && itera >= MIN_ITERATIONS) {
      print(itera)
      print(delta)
      print(max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera]))
    }
    if(delta < MIN_CHANGE || (itera >= MIN_ITERATIONS &&
                              max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera]) < min(epsilon_2))) {
      break
    }
    if(F) { # SCRATCH TEMP
      precision <- abs(solve(D_mat) %*% d_vec)
      if(itera >= MIN_ITERATIONS)
        print(paste(precision[28],max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera])))
    }
  }
  #plot(log(pred_decrease_vals, 10))
  #lambda[lambda < lam_threshold] <- lam_threshold
  #d1_d2 <- get_d1_d2(Y, lambda)
  H <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2)$H
  eps = 1e-6
  if(min(eigen(H)$values) < eps) {
    I <- solve(psd(H, epsilon = eps))
  } else {
    I <- solve(H)
  }
  precision <- abs(solve(D_mat) %*% d_vec)
  precision[is.na(diag(I))] <- PRECISION.THRESHOLD + 100
  converged_vec <- check_converged_vec(X1,X2,my_beta, itera, n.iter,
                                       error_vec, precision, PRECISION.THRESHOLD)
  names(error_vec) <- colnames(my_beta)
  return(list(alpha1 = alpha1, alpha2 = alpha2, converged = any(converged_vec), I = I, d = d_vec_o,
              n.iter = itera, log_l = prev_ll, precision = precision, prediction = lambda,
              converged_vec = converged_vec, error_vec = error_vec))
}

# KEY FUNCTION 1
get_de_gene_fits <- function(X1,X2,my_beta, nUMI, gene_list, cell_types, puck, barcodes, sigma_init, test_mode,
                             numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.05, params_to_test = 2, logs=F) {
  results_list <- fit_de_genes(X1,X2,my_beta, nUMI, gene_list, puck, barcodes,
                               sigma_init, test_mode, numCores = numCores,
                               sigma_gene = sigma_gene,
                               PRECISION.THRESHOLD = PRECISION.THRESHOLD, logs = logs)
  N_genes <- length(results_list)
  intercept_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  mean_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  all_vals <- array(0, dim = c(N_genes, dim(X2)[2],length(cell_types)))
  dimnames(all_vals)[[1]] <- gene_list
  dimnames(all_vals)[[3]] <- cell_types
  con_val <- logical(N_genes)
  ll_val <- numeric(N_genes)
  n_val <- numeric(N_genes)
  sigma_g <- numeric(N_genes)
  names(sigma_g) <- gene_list
  I_val <- list()
  names(n_val) <- gene_list
  names(con_val) <- gene_list
  names(ll_val) <- gene_list
  rownames(mean_val) <- gene_list; colnames(mean_val) <- cell_types
  rownames(intercept_val) <- gene_list; colnames(intercept_val) <- cell_types
  d_vals <- matrix(0,nrow=N_genes,ncol=dim(X2)[2]*length(cell_types))
  s_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  precision_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_all <- matrix(FALSE, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  error_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  rownames(precision_mat) <- gene_list; rownames(con_all) <- gene_list
  rownames(s_mat) <- gene_list; rownames(con_mat) <- gene_list; rownames(error_mat) <- gene_list
  colnames(s_mat) <- get_param_names(X1,X2, cell_types)
  colnames(precision_mat) <- get_param_names(X1,X2, cell_types)
  colnames(con_all) <- get_param_names(X1,X2, cell_types)
  colnames(con_mat) <- cell_types
  colnames(error_mat) <- cell_types
  rownames(d_vals) <- gene_list
  for(i in 1:N_genes) {
    sigma_g[i] <- results_list[[i]]$sigma_s_best
    res <- results_list[[i]]$res
    d_vals[i,] <- res$d
    mean_val[i,] <- res$alpha2[params_to_test[1],]
    intercept_val[i,] <- res$alpha2[1,]
    all_vals[i, ,] <- res$alpha2
    con_val[i] <- res$converged
    precision_mat[i,] <- res$precision
    ll_val[i] <- res$log_l
    n_val[i] <- res$n.iter
    I_val[[i]] <- res$I
    s_mat[i,] <- sqrt(diag(I_val[[i]]))
    con_mat[i,] <- res$converged_vec
    con_all[i,] <- res$precision < PRECISION.THRESHOLD
    error_mat[i,] <- res$error_vec
  }
  return(list(mean_val = mean_val, con_val = con_val, ll_val = ll_val, I_val = I_val, s_mat = s_mat,
              n.iter = n_val,d_vals = d_vals, intercept_val = intercept_val, all_vals = all_vals,
              precision_mat = precision_mat, sigma_g = sigma_g, con_mat = con_mat, con_all = con_all, error_mat = error_mat))
}

estimate_gene_wrapper <- function(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = T, PRECISION.THRESHOLD = 0.05) {
  if(sigma_gene)
    return(choose_sigma_gene(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD))
  else {
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    return(list(sigma_s_best = -1, res = res))
  }
}

estimate_effects_trust <- function(Y, X1, X2, my_beta, nUMI, test_mode, verbose = F,
                                   n.iter = 200, MIN_CHANGE = 1e-3, PRECISION.THRESHOLD = 0.05,
                                   alpha1_init = NULL, alpha2_init = NULL, MIN_ITERATIONS = 15) {
  my_beta<- sweep(my_beta,1, nUMI, '*')
  solveIRWLS.effects_trust(Y,X1,X2, my_beta, test_mode, verbose = verbose,
                           n.iter = n.iter, MIN_CHANGE = MIN_CHANGE,
                           PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                           alpha1_init = alpha1_init, alpha2_init = alpha2_init,
                           MIN_ITERATIONS = MIN_ITERATIONS)
}