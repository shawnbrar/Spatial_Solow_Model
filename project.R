library(data.table)
library(Matrix)
library(spdep)
library(spatialreg)

source("functions.R") # load functions to be used
## reading GDP Data
gdp <- fread("data/gdpc.csv", select = c("geo", "TIME_PERIOD", "OBS_VALUE"))
## reading Capital Data
cap <- fread("data/capc.csv", select = c("geo", "TIME_PERIOD", "OBS_VALUE"))
## reading Population Data
pop <- fread("data/population", select = c("age", "geo", "TIME_PERIOD", "OBS_VALUE"))
pop <- pop[, sum(OBS_VALUE), by = c("geo", "TIME_PERIOD")]
## loading HDI Data
load("data/hdi.RData") ## prepared in hdi_prepare.R
# Load gerd data
gerd_data <- fread("data/rd_e_gerdtot__custom_3182159_linear.csv", select = c("geo", "TIME_PERIOD", "OBS_VALUE"))
gerd_data2 <- gerd_data[!duplicated(gerd_data[, c("geo", "TIME_PERIOD")]), ]
# Combing capital and population datasets into gdp
gdp <- gdp[cap, on = c(geo = "geo", TIME_PERIOD = "TIME_PERIOD"), nomatch = 0]
gdp <- gdp[pop, on = c(geo = "geo", TIME_PERIOD = "TIME_PERIOD"), nomatch = 0]
gdp <- gdp[gerd_data2, on = c(geo = "geo", TIME_PERIOD = "TIME_PERIOD"), nomatch = 0]
colnames(gdp) <- c("geo", "year", "gdp", "cap", "population", "gerd")
## combining gdp with hdi
cc <- fread("data/country_codes_V202201.csv")
gdp[cc, c("country", "iso3_code") := .(i.country_name_abbreviation, i.iso_3digit_alpha), on = c(geo = "iso_2digit_alpha")]
gdp <- gdp[hdi2[, c(1, 3, 4)], on = c(iso3_code = "iso3", year = "year"), nomatch = 0]
rm(cap, pop, hdi2)
gdp[, population := population/1000000]
gdp[, c("gdp", "cap") := .(gdp/population, cap/population)] # gdp and cap per
                                                            # labour employed
# Create time lags for gdp and population
gdp[, c("gdp_t1", "pop_t1") := lapply(.SD, shift), by = "geo", .SDcols = c("gdp", "population")]
# Calculate population growth rate and add 0.05
gdp[, pop_gr := ((population - pop_t1)/pop_t1) + 0.05]
gdp[, pop_t1 := NULL]
# Calculate share of investment in gdp i.e. savings rate
gdp[, sav := cap/gdp]
gdp <- na.omit(gdp)
# Create logs of the variables
cols <- colnames(gdp)[c(3:4, 6, 9:10, 12)]
gdp[, (cols) := lapply(.SD, log), .SDcols = cols]
gdp[, gdp_diff := gdp - gdp_t1]
#lapply(1:11, function(x) any(is.infinite(gdp[[x]]))) # check if any value in the data is infinite
## Make the spatial weights matrix ####
cc <- fread("data/country_codes_V202201.csv")
capitals <- rgdal::readOGR("./data/shp/", "eurcap", verbose = FALSE)
eu_countries <- c("Austria", "France", "Malta", "Belgium", "Germany",
                  "Netherlands", "Bulgaria", "Poland", "Croatia",
                  "Hungary", "Portugal", "Cyprus", "Ireland", "Romania",
                  "Czechia", "Italy", "Slovakia", "Denmark", "Latvia",
                  "Slovenia", "Estonia", "Lithuania", "Spain", "Finland",
                  "Luxembourg", "Sweden")
capitals@data[capitals@data$CNTRY_NAME == "Czech Republic" , 3] <- "Czechia"
capitals <- capitals[grep("National", capitals@data$STATUS), ]
capitals <- capitals[capitals@data$CNTRY_NAME %in% eu_countries, ]
coords <- as.data.table(coordinates(capitals))
coords <- cbind(coords, capitals@data[, c(1, 3)])
coords[cc[, c(2, 4, 5)], c("iso2", "iso3") := .(i.iso_2digit_alpha, i.iso_3digit_alpha), on = c(CNTRY_NAME="country_name_abbreviation")]
gdp <- gdp[geo %in% coords$iso2, ]
#gdp <- gdp[year %in% 1981:2021]
# construct distance matrix
years <- sort(unique(gdp$year))
matrices <- lapply(years, function(x) 1/as.matrix(dist(coords[iso2 %in% gdp[year == x, ]$geo, 1:2])))
spt_temp_mat <- bdiag(matrices)
diag(spt_temp_mat) <- 0L
spt_temp_mat <- spt_temp_mat/rowSums(spt_temp_mat)
means <- apply(spt_temp_mat, 1, mean)
spt2 <- (spt_temp_mat >= rowMeans(
  replace(spt_temp_mat,spt_temp_mat == 0, NA), na.rm = TRUE))*1 # Distance Matrix
spt2 <- spt2/rowSums(spt2)
rm(matrices, spt_temp_mat, years)
gc()

# Construct Bilateral Trade Flow Matrix
load("data/export_imports.RData") # prepared in imp_exp_data.R
exp_all[cc[, c(1, 2, 4)],
        c("exp_country_name", "exp_geo") := .(i.country_name_abbreviation, i.iso_2digit_alpha),
        on = c(i="country_code")]
exp_all[cc[, c(1, 2, 4)],
        c("imp_country_name", "imp_geo") := .(i.country_name_abbreviation, i.iso_2digit_alpha),
        on = c(j="country_code")]
setorder(exp_all, "exp_geo", "imp_geo")
years <- unique(exp_all$year)
gdp_exp <- copy(gdp)
gdp_exp <- split(gdp_exp, gdp_exp$year) ## split the datatable by year, so we can
                                        ## reduce the datatable for each year depending
                                        ## on each year imports and exports
gdp_exp[!(names(gdp_exp) %in% as.character(1996:2019))] <- NULL # put 1981 till 1995 as null

## we will have 1996:2021 in gdp and 1995:2020 in exp_all
years <- as.integer(names(gdp_exp))
exp_all <- exp_all[year %in% (years-1), ]
system.time(mat_data <- lapply(1:length(years), function(x) exp_mat(years[x], exp_all, gdp_exp[[x]])))
exp_matrix <- bdiag(sapply(mat_data, "[[", 1)) ## export matrix
exp_matrix <- exp_matrix / rowSums(exp_matrix)
imp_matrix <- t(exp_matrix)
imp_matrix <- imp_matrix / rowSums(imp_matrix)
bi_trd_fl_mat <- imp_matrix + exp_matrix
gdp_exp <- rbindlist(lapply(mat_data, "[[", 2)) ## export data
# creating geodist matrix to be multiplied to bilateral trade flow matrices
years <- sort(unique(gdp_exp$year))
matrices <- lapply(years, function(x) 1/as.matrix(dist(coords[iso2 %in% gdp_exp[year == x, ]$geo, 1:2])))
spt_temp_mat <- bdiag(matrices)
diag(spt_temp_mat) <- 0L
spt_temp_mat <- spt_temp_mat/rowSums(spt_temp_mat)
means <- apply(spt_temp_mat, 1, mean)
spt3 <- (spt_temp_mat >= rowMeans(
  replace(spt_temp_mat,spt_temp_mat == 0, NA), na.rm = TRUE))*1 # Distance Matrix
geodist_export <- spt3 * exp_matrix
geodist_import <- spt3 * imp_matrix
geodist_bilat <- spt3 * bi_trd_fl_mat
rm(cc, exp_all, mat_data, years)

## create a list of weight matrices
spt_wgt_mat_list <- list("spt2" = spt2, "imp_mat" = imp_matrix, "exp_mat" = exp_matrix, "bi_trd_fl" = bi_trd_fl_mat,
                         "geodist_import" = geodist_import, "geodist_export" = geodist_export, "geodist_bilat" = geodist_bilat)
listw_obj <- lapply(spt_wgt_mat_list, mat2listw, style = "W")
rm(spt2, imp_matrix, exp_matrix, bi_trd_fl_mat, geodist_bilat, geodist_export,
   geodist_import);gc()

## Create WX ####
wX <- as.data.table(as.matrix(spt_wgt_mat_list[[1]] %*% as.matrix(gdp[, c(6, 9, 11:12)])))
gdp[, c("w_gerd", "w_hdi", "w_pop_gr", "w_sav") := .(wX$gerd, wX$hdi, wX$pop_gr, wX$sav)]
wX <- as.data.table(as.matrix(spt_wgt_mat_list[["bi_trd_fl"]] %*% as.matrix(gdp_exp[, c(6, 9, 11:12)]))) ## bilateral trade
gdp_exp[, c("w_gerd", "w_hdi", "w_pop_gr", "w_sav") := .(wX$gerd, wX$hdi, wX$pop_gr, wX$sav)]

## Start the modelling process ####
# given a weight matrix calculate the four models - idea #
coefs_to_find <- c("(Intercept)", "gdp_t1",  "sav",  "hdi", "pop_gr", "gerd", "w_sav",  "w_hdi", "w_pop_gr", "w_gerd", "lag.sav", "lag.hdi", "lag.gerd", "lag.pop_gr")

geodist_models <- model_obj(listw_obj[[1]], gdp)
import_models <- model_obj(listw_obj[[2]], gdp_exp)
export_models <- model_obj(listw_obj[[3]], gdp_exp)
bi_trd_fl_models <- model_obj(listw_obj[[4]], gdp_exp)
geodist_import_models <- model_obj(listw_obj[[5]], gdp_exp)
geodist_export_models <- model_obj(listw_obj[[6]], gdp_exp)
geodist_bitrd_models <- model_obj(listw_obj[[7]], gdp_exp)

all_models_obj <- c(geodist_models, import_models, export_models, bi_trd_fl_models,
                    geodist_import_models, geodist_export_models, geodist_bitrd_models)
names(all_models_obj) <- paste0(rep(c("geodist_", "import_", "export_", "bi_trd_fl_",
                                      "geodist_import_", "geodist_export_", "geodist_bitrd_"), each = 4),
                                names(all_models_obj))
summ_model_obj <- lapply(all_models_obj, summary)
rm(geodist_models, import_models, export_models, bi_trd_fl_models,
   geodist_bitrd_models, geodist_export_models, geodist_import_models);gc()
ext_mod_coefs <- lapply(1:length(summ_model_obj),
                        function(x) extract_coefs(summ_model_obj[[x]],
                                                  all_models_obj[[x]])) # extracted model coefficients
names(ext_mod_coefs) <- names(all_models_obj)
rm(summ_model_obj);gc()



LR.Sarlm(all_models_obj[[3]], all_models_obj[[2]])
#LR Tests - SDM vs SAR - SAC vs SEM - SAC vs SAR # data -
#all_models_obj[grep("geodist", names(all_models_obj))]
spatialreg::test


## LM Tests ####
suppressWarnings({lmtests_obj <- lapply(c("geodist", "import", "export", "bi_trd_fl", "geodist_import", "geodist_export", "geodist_bitrd"),
       function(x) lmtests(all_models_obj[grep(x, names(all_models_obj))], x))})
names(lmtests_obj) <- c("geodist", "import", "export", "bi_trd_fl", "geodist_import", "geodist_export", "geodist_bitrd")

lmtest_res_obj <- lapply(c("geodist", "import", "export", "bi_trd_fl", "geodist_import", "geodist_export", "geodist_bitrd"),
       function(x) lmrestlt(lmtests_obj[[x]], x))
