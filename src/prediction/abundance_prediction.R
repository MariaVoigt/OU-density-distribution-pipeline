#-----------------------------------------------#
# Model prediction:                             #
# the next part is for making model predictions #
# throughout study area                         #
# define predictors of 'grid'                   #
#-----------------------------------------------#

#----------------#
# Load Libraries #
#----------------#

suppressMessages(library(parallel))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(rgdal))
suppressMessages(library(sp))
suppressMessages(library(rgeos))
suppressMessages(library(raster))
suppressPackageStartupMessages(library(optparse))

cl <- makeForkCluster(outfile = "")
registerDoParallel(cl)

#------------------------#
# command line arguments #
#------------------------#
print(paste(Sys.time(), "0. start run"))
args <- commandArgs(trailingOnly = TRUE)
print(paste("args", args))
indir <- args[1]
print(paste("indir ", indir))
outdir <- args[2]
print(paste("outdir ", outdir))
indir_predictors <- args[3]
print(paste("indir predictors", indir_predictors))
year_to_predict <- as.numeric(args[4])
print(paste("year " , year_to_predict))

indir_fun <- "~/orangutan_density_distribution/src/functions/"
crs_aea <- "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0
  +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # final projection


#-------------------------------#
# Load and prepare coefficients #
#-------------------------------#
# Load functions
source(file.path(indir_fun, "generic/path.to.current.R"))
print("function loaded")

# Load coefficients and weights
abundMod_results_path <- path.to.current(indir, "abundMod_results", "rds" )
print(paste("this is abundMod_results_path:", abundMod_results_path))
abundMod_results <- readRDS(abundMod_results_path)
# exclude the first column, which contains models, exclude the coefficient of the
# autocorellation term and the weighted aic of the model
# we are excluding the autocorellation term, because the mean is zero, and thus
# it will not have an influence on the overall outcome
coeffs <- abundMod_results[ ,  grepl(x=colnames(abundMod_results),
                                     pattern="coeff")]
# all predictors that are not included in a specific model have NA as their
# coefficient, which is replaced by 0, so that the estimate is also 0
coeffs[is.na(coeffs) == T] <- 0

# Load estimates for observation and grid
# these are the predictors that will be used in the prediction
# THEY MUST BE IN THE ORDER IN WHICH THEY APPEAR IN THE COEFFICIENTS
# CODE THIS
# here only separate after first _
predictor_names_coeffs <- gsub("coeff_","", names(coeffs))
#UNDERSTAND HERE WHAT IS HAPPENING
#interaction_terms_names <- predictor_names_coeffs[predictor_names_coeffs %in%
#                           paste0("year:", predictor_names_coeffs)]
#interaction_terms_names <- gsub("year:", "", interaction_terms_names)
#quadratic_terms_names <- predictor_names_coeffs[predictor_names_coeffs %in%
#                                                 paste0("I(", predictor_names_coeffs, "^2)")]
#quadratic_terms_names <- gsub("I\\(|\\^2\\)", "", quadratic_terms_names )

quadratic_terms_names <- c("rain_dry")
predictor_names_coeffs <- predictor_names_coeffs[predictor_names_coeffs != "(Intercept)"]
# don't include interaction or quadratic term
predictor_names <- predictor_names_coeffs[!grepl("I(*)", predictor_names_coeffs)]

#----------------------------#
# Load and prepare estimates #
#----------------------------#
print(paste("these are predictor names: ", predictor_names))
predictors_path <- path.to.current(indir_predictors, paste0("predictors_abundance_",
                                                 year_to_predict),"rds")
print(paste("this is predictors path", predictors_path))
predictors <- readRDS(predictors_path) %>%
 dplyr::filter(predictor %in% predictor_names) %>%
    dcast(id + z_year ~ predictor,  value.var = "scaled_value")

predictors$year <- predictors$z_year
predictors$z_year <- NULL
str(predictors)



print(paste("this is nrow predictors", nrow(predictors)))
#--------------------------#
# PREDICTION FOR each year #
#--------------------------#
#predictor_estimates must comprise covariates and dumy codes in the sequence in which they appear
# in the coefficients
# HERE PAY ATTENTION FOR INTERACTIONS
# AND IF THERE IS MORE THAN ONE QUADRATIC TERM
intercept <- rep(1, nrow(predictors))
predictor_estimates <- cbind( intercept,
                    predictors[ , predictor_names],
                    predictors[ ,quadratic_terms_names] * predictors[ ,quadratic_terms_names])

names(predictor_estimates) <- c("intercept", predictor_names,
                                paste0("I(", quadratic_terms_names, "^2)"))



# Alternative here to loop through id, which is added to predictor_estimates
# more correct in terms of being sure that the right thing is done,
# but takes a bit longer (52s, to 35s for 100 rows)
## PLUS PAY ATTENTION, IF PREDICTIONS NOT SAME NROW--> VALUES GET RECYCLED

print(paste("1. start pred_per_cell", Sys.time()))
pred_per_cell <- foreach(i = 1:nrow(predictor_estimates), .combine = c)  %dopar% {
# pred_per_cell <- foreach(i = 1:100, .combine = c)  %dopar% {
t_predictor_estimates <- t( predictor_estimates[i, ])
pred_estimates_wcoeffs  <- data.frame(mapply(`*`, coeffs, t_predictor_estimates, SIMPLIFY = F))
pred_estimates_sum <- apply(pred_estimates_wcoeffs, 1, sum)
    pred_estimates_weighted <- pred_estimates_sum * abundMod_results$w_aic
    pred_estimates_calc <- sum(pred_estimates_weighted)
return(exp(pred_estimates_calc))
}

print(paste(Sys.time(), "2. finished dopar loop"))

# is this correct -> ????
pred_per_cell <- as.data.frame(cbind(predictors$id, pred_per_cell))
names(pred_per_cell) <- c("id", "abundance_pred")

# exclude NAN
pred_per_cell <- pred_per_cell[!is.nan(pred_per_cell$abundance_pred), ]


print(paste(Sys.time(), "sum predicted for ", year_to_predict,
            sum(pred_per_cell$abundance_pred)))
print(paste(Sys.time(), "range predicted for ", year_to_predict,
            range(pred_per_cell$abundance_pred)))
print(paste(Sys.time(), "density predicted for ", year_to_predict,
            mean(pred_per_cell$abundance_pred)))
print(paste(Sys.time(), "nr of values over 10 ", year_to_predict,
            sum(pred_per_cell$abundance_pred > 10)))

saveRDS(pred_per_cell,
        file = file.path(outdir,
                         paste0("abundance_pred_per_cell_",
                                year_to_predict,"_",
                                Sys.Date(), ".rds")))

save.image(file.path(outdir, "image_before_map.RData"))

#-----------------------#
# convert output to map #
#-----------------------#
print(paste(Sys.time(), "3. Start making map"))

geography_grid_path <- path.to.current(indir_predictors,
                                       paste0("geography_", year_to_predict), "rds")
geography_grid <- readRDS(geography_grid_path)

geography_grid_for_join <- dplyr::select(geography_grid, id, x_start, y_start)

pred_per_cell_sp <- left_join(geography_grid_for_join,
                                  pred_per_cell, by = "id") %>%
                        as.data.frame()


pred_per_cell_sp <- SpatialPointsDataFrame(coords =
                               cbind(pred_per_cell_sp$x_start,
                                     pred_per_cell_sp$y_start),
                                  data = pred_per_cell_sp,
                                  proj4string = CRS(crs_aea),  match.ID = T)

writeOGR(pred_per_cell_sp, dsn = outdir,
         layer = paste0("prediction_shp_", year_to_predict), driver = "ESRI Shapefile")
# pay attention that the tif gets a new name because later we delete all
# with the same name as the shapefile
system_string <- paste0("gdal_rasterize -ot 'Float32' -a abndnc_ ",
                        "-tr 1000.0 1000.0 -l prediction_shp_",
                         year_to_predict," ", outdir, "/",
                        "prediction_shp_", year_to_predict, ".shp",
                        " ", outdir, "/",
                        "prediction_map_", year_to_predict, ".tif")
system(system_string)
#delete the shapefile again because clogs up my system
system_string_del <- paste0("rm -f ", outdir, "/prediction_shp_",
                            year_to_predict, ".*")
system(system_string_del)



save.image(file.path(outdir, paste0("abundance_pred_image_", year_to_predict, "_",
                                    Sys.Date(), ".RData")))



print(paste(Sys.time(), "4. wrote results and done :-)"))
