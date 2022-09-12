## Functions for the project

exp_mat <- function(yearz, exp_data, gdp_data){
  # function to make the export matrix and export data
  dexp <- exp_data[year == (yearz - 1), ][exp_geo %in% gdp_data$geo, ][imp_geo %in% gdp_data$geo, ]
  dexp[, c("exp_rank", "imp_rank") := .(frank(exp_geo, ties.method = "dense"), frank(imp_geo, ties.method = "dense"))]
  mat <- sparseMatrix(i = dexp$exp_rank, j = dexp$imp_rank, x = dexp$value)
  not_in <- unique(gdp_data$geo)[which(!(unique(gdp_data$geo) %in% unique(dexp$exp_geo)))]
  return(list(mat, gdp_data[!(geo %in% not_in), ]))
}

model_obj <- function(listw_use, data){
  # function to create the four models
  sar <- lagsarlm(gdp_diff ~ gdp_t1 + sav + hdi + gerd + pop_gr + as.factor(geo), data, listw_use)
  sem <- errorsarlm(gdp_diff ~ gdp_t1 + sav + hdi + gerd + pop_gr + as.factor(geo), data, listw_use, etype = "error")
  sdm <- lagsarlm(gdp_diff ~ gdp_t1 + sav + hdi + gerd + pop_gr + as.factor(geo), data, listw_use, Durbin = ~sav + hdi + gerd + pop_gr) #+ w_sav + w_hdi + w_gerd + w_pop_gr
  slx <- lmSLX(gdp_diff ~ gdp_t1 + sav + hdi + gerd + pop_gr + as.factor(geo), data, listw_use, Durbin = ~sav + hdi + gerd + pop_gr)
  return(list("sar_model" = sar, "sem_model" = sem, "sdm_model" = sdm, "slx_model" = slx))
}

extract_coefs <- function(model_sum, model_obj){
  # function to extract model coefficients
  a <- rownames(model_sum$coefficients)
  b <- names(model_sum$coefficients)
  coefs <- data.table(coefs = if(is.null(a)) {b} else {a}, coef(model_sum))
  coefs <- coefs[coefs %in% coefs_to_find, ]
  vec <- c()
  if(is.null(model_sum[["rho"]]) == FALSE){
    rho <- list("rho", model_sum[["rho"]], model_sum[["rho.se"]], model_sum[["rho"]]/model_sum[["rho.se"]], 2 * (1 - pnorm(abs(model_sum[["rho"]]/model_sum[["rho.se"]]))))
    vec[length(vec)+1] <- "Yes"
  } else {
    vec[length(vec)+1] <- "No"
  }
  if(is.null(model_sum[["lambda"]]) == FALSE){
    lambda <- list("lambda", model_sum[["lambda"]], model_sum[["lambda.se"]], model_sum[["lambda"]]/model_sum[["lambda.se"]], 2 * (1 - pnorm(abs(model_sum[["lambda"]]/model_sum[["lambda.se"]]))))
    vec[length(vec)+1] <- "Yes"
  } else {
    vec[length(vec)+1] <- "No"
  }
  loglikeli <- list("Log Lik", as.numeric(logLik(model_obj)), 1, 1, 1)
  aic <- list("AIC", AIC(model_obj), 1, 1, 1)
  coefs <- rbindlist(list(coefs, if(vec[1] == "Yes") rho, if(vec[2] == "Yes") lambda, loglikeli, aic))
  return(coefs)
}

lmtests <- function(dataz, matrix_type){
  # function to perform lmtests
  q <- paste0(matrix_type, c("_sar_", "_sdm_"), "model")
  r <- paste0(matrix_type, c("_sem_", "_sdm_"), "model")
  t <- paste0(matrix_type, c("_slx_", "_sdm_"), "model")
  qsd <- list(q, r, t)
  lis <- lapply(qsd, function(x) lmtest::lrtest(dataz[[x[1]]], dataz[[x[2]]]))
  return(lis)
}

lmrestlt <- function(list_obj, mat_type){
  # Function to crerate lm test result for 3 tests
  retz <- data.table(Mat = rep(mat_type, times = 3), Res = c("SAR", "SEM", "SLX"),
                     UnRes = c("SDM", "SDM", "SDM"), chi_stat = as.numeric(sapply(list_obj, "[[", 4)[2,]),
                     prob = as.numeric(sapply(list_obj, "[[", 5)[2,]))
  return(retz)
}
