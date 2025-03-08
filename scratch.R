












tslope_effs_mask <- grepl("tslope", merfish.laminar.model$param.names) & grepl("beta", merfish.laminar.model$param.names)
bs_fitted_params <- merfish.laminar.model$bs_params
tslope_effects <- c(bs_fitted_params[,tslope_effs_mask])
hist(tslope_effects)
sd(tslope_effects)
mean(tslope_effects)

mask <- grepl("beta_tslope", testnames)
plot(test[mask])
max(test[mask])



pools <- c()
for (i in 1:length(merfish.laminar.model[["token.pool"]])) {
  l <- length(merfish.laminar.model[["token.pool"]][[i]])
  if (l > 0) {
    pools <- c(pools, l)
  }
}
hist(pools)
mean(pools)
counts <- merfish.laminar.model$count.data.summed
resamples <- merfish.laminar.model$resample.demo

mask <- counts$child == "Rorb" & counts$ran == "none" & counts$treatment == "ref"
plot(
  counts$bin[mask], 
  counts$count[mask], 
  pch = 19, 
  col = "blue", 
  xlab = "Bin", 
  ylab = "Count")
j <- sample(1:ncol(resamples), 1)
points(
  counts$bin[mask], 
  resamples[mask, j], 
  pch = 19, 
  col = "red")








pn <- merfish.laminar.model$param.names
ps <- merfish.laminar.model$fitted.parameters
ps <- merfish.laminar.model$fitted.parameters
mc <- sum(counts$count.log[counts$ran != "none"], na.rm = TRUE)/4
mean_wf <- c()
sd_wf <- c()
mean_wf_obs <- c()
for (r in 1:4) {
  param_mask <- grepl(paste0("wfactor_rate_", r, "_X"), pn)
  hist(ps[param_mask])
  print(mean(ps[param_mask]))
  print(sum(counts$count.log[counts$ran == as.character(r)], na.rm = TRUE)/mc - 1)
  print("------")
  mean_wf <- c(mean_wf, mean(ps[param_mask]))
  sd_wf <- c(sd_wf, sd(ps[param_mask]))
  mean_wf_obs <- c(mean_wf_obs, sum(counts$count.log[counts$ran == as.character(r)], na.rm = TRUE)/mc - 1)
}
ds <- c()
ds_max <- c()
for (r in 1:4) {
  
  d <- dnorm(mean_wf_obs[r], mean_wf[r], sd_wf[r]/sqrt(4))
  ds <- c(ds, d)
  d_max <- dnorm(mean_wf[r], mean_wf[r], sd_wf[r]/sqrt(4))
  ds_max <- c(ds_max, d_max)
  
}
results <- data.frame(mean_wf_obs, mean_wf, sd_wf, ds, ds_max)



plot.decomposition(
  merfish.laminar.model,
  child = "Rorb"
)






