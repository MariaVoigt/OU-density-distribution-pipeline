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
suppressMessages(registerDoParallel(cl))


#-----------------------------#
# command line option parsing #
#-----------------------------#
print(paste("Start abundance_prediction script", Sys.time()))


option_list <- list (
  make_option(c("-i", "--input-directory"),  dest = "input_directory",
              type = "character", help = "directory with input files",
              metavar = "/path/to/input-dir"),
  make_option(c("-o", "--output-directory"), dest = "output_directory",
               type = "character", help = "directory with output files",
               metavar = "/path/to/output-dir"),
  make_option("--exclude-year",    dest = "exclude_year", type = "integer",
              default = NA, help = "year to exclude", metavar = "2015"),
  make_option("--year-to-predict",
              dest = "year_to_predict",
              type = "integer",
              help = "year of the survey years (1994:2015) to predict abundance to",
              metavar = "2015"),
  make_option("--focal-change-predictor",
              dest = "focal_change_predictor",
              type = "character",
              help = "predictor that is allowed to vary in stability analysis",
              metavar = "year"),
  make_option(c("-q", "--quiet"), dest = "verbose_script",
              action = "store_false",
              default = TRUE,
              help = "don't print all intermediate results"))

# verbose option a bit counterintuitive
# because I make action store_false, when I say -q that
# means that verbose == F, which is quiet


options <- parse_args(OptionParser(option_list=option_list))

# sanity checks

if (is.null(options$input_directory)) {
  stop("input directory not specified, check --help")
}


if (is.null(options$output_directory)) {
  stop("output directory not specified, check --help")
}

exclude_year_possibilities <- c(1999:2015)

if (!is.na(options$exclude_year) && !(options$exclude_year %in% exclude_year_possibilities)) {
  stop(paste("exclude year must be between", min(exclude_year_possibilities),
             "and", max(exclude_year_possibilities)))
}

year_to_predict_possibilities <- c(1999:2015)

if (!(options$year_to_predict %in% year_to_predict_possibilities)) {
  stop(paste("the year to predict to year must be between",
             min(year_to_predict_possibilities),
             "and", max(year_to_predict_possibilities)))
}


# is quiet?
is_verbose <- options$verbose_script
# input directory
indir <- options$input_directory
if(is_verbose){print(paste("indir", indir))}


# directory in which output is written
outdir <- options$output_directory
if(is_verbose){print(paste("outdir", outdir))}

year_to_predict <- as.numeric(options$year_to_predict )
print(paste("for year to predict" , year_to_predict))


exclude_year <- as.numeric(options$exclude_year)
if(is_verbose){print(paste("exclude year", exclude_year))}

if(is_verbose){print(paste(Sys.time(), "0. start run"))}
#------------------------#
# command line arguments #
#------------------------#

indir_fun <- "~/orangutan_density_distribution/src/functions/"

crs_aea <- "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0
  +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # final projection

if (is.na(exclude_year)){
  name_suffix <- ""} else {
    name_suffix <- paste0(exclude_year, "_")
  }

#----------------#
# Load functions #
#----------------#
# Load functions
source(file.path(indir_fun, "generic/path.to.current.R"))
source(file.path(indir_fun, "project_functions/scale.predictors.R"))


#----------------------#
# Prepare coefficients #
#----------------------#
# Load coefficients and weights
abundMod_results_path <- path.to.current(indir, "abundMod_results", "rds" )
if(is_verbose){print(paste("this is abundMod_results_path:", abundMod_results_path))}
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
geography_path <- path.to.current(indir, paste0("geography_",
                                                            year_to_predict),"rds")
if(is_verbose){print(paste("this is geography path", geography_path))}
geography <- readRDS(geography_path)

predictors_path <- path.to.current(indir, paste0("predictors_abundance_",
                                                 year_to_predict),"rds")
if(is_verbose){print(paste("this is predictors path", predictors_path))}
predictors <- readRDS(predictors_path)

# import here the already z-transformed or not
predictors_obs_path <- path.to.current(indir,
                                       paste0("predictors_observation_scaled_",
                                              name_suffix),
                                              "rds")
if(is_verbose){print(paste("this is predictors-obs path", predictors_obs_path))}

predictors_obs <- readRDS(predictors_obs_path)

# Scale the grid-predictors using mean and sd of predictors_obs

# predictors for scaling
predictor_names_for_scaling <- c( "dem", "slope", "temp_mean", "rain_dry", "rain_var",
                                  "ou_killings", "ou_killing_prediction", "human_pop_dens",
                                  "perc_muslim", "peatswamp", "lowland_forest", "lower_montane_forest" ,
                                  "road_dens", "distance_PA", "fire_dens", "deforestation_hansen",
                                  "deforestation_gaveau", "plantation_distance", "pulp_distance", "palm_distance",
                                  "dom_T_OC", "dom_T_PH")

predictor_names_add <- c("year", "x_center", "y_center")

predictors <- rename(predictors, unscaled_value = value)

# calculate x and y center
geography$unscaled_x_center <-  rowMeans(cbind(geography$x_start, geography$x_end), na.rm = T)
geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = T)



# function here
predictors_grid <- scale.predictors.grid(predictor_names_for_scaling,
                                  predictor_names_add,
                                  predictors,
                                  predictors_obs,
                                  geography)


saveRDS(predictors_grid, file.path(outdir, paste0("predictors_grid_scaled_", name_suffix,
                                 year_to_predict, "_", Sys.Date(), ".rds")))

#-------------------------------------#
# Prepare analysis of decline drivers #
#-------------------------------------#
if(!is.na(focal_change_predictor)){
  # load 1999 value to add the values, but use the scaled table
  predictors_grid_1999_path <- path.to.current(indir,
                                         paste0("predictors_grid_scaled_",
                                                name_suffix, 1999),
                                         "rds")
  predictors_grid_1999 <- readRDS(  predictors_grid_1999_path)
  change_predictors <- c("year", "peatswamp",
                         "lowland_forest",
                         "lower_montane_forest",
                         "deforestation_hansen")
  change_predictors <- change_predictors[!change_predictors %in% focal_change_predictor]
  for (change_predictor in change_predictors){
    # go through all predictors and set their value to the 1999 value
    predictors_grid[ , change_predictor] <- predictors_grid_1999[ , change_predictor]
  }
}

#--------------------------#
# PREDICTION FOR each year #
#--------------------------#
#predictor_estimates must comprise covariates and dumy codes in the sequence in which they appear
# in the coefficients
# HERE PAY ATTENTION FOR INTERACTIONS
# AND IF THERE IS MORE THAN ONE QUADRATIC TERM
intercept <- rep(1, nrow(predictors_grid))
predictor_estimates <- cbind( intercept,
                              predictors_grid[ , predictor_names],
                              predictors_grid[ ,quadratic_terms_names] * predictors_grid[ ,quadratic_terms_names])

names(predictor_estimates) <- c("intercept", predictor_names,
                                paste0("I(", quadratic_terms_names, "^2)"))



# Alternative here to loop through id, which is added to predictor_estimates
# more correct in terms of being sure that the right thing is done,
# but takes a bit longer (52s, to 35s for 100 rows)
## PLUS PAY ATTENTION, IF PREDICTIONS NOT SAME NROW--> VALUES GET RECYCLED

if(is_verbose){print(paste("1. start pred_per_cell", Sys.time()))}
pred_per_cell <- foreach(i = 1:nrow(predictor_estimates), .combine = c)  %dopar% {
  # pred_per_cell <- foreach(i = 1:100, .combine = c)  %dopar% {
  t_predictor_estimates <- t( predictor_estimates[i, ])
  pred_estimates_wcoeffs  <- data.frame(mapply(`*`, coeffs, t_predictor_estimates, SIMPLIFY = F))
  pred_estimates_sum <- apply(pred_estimates_wcoeffs, 1, sum)
  pred_estimates_weighted <- pred_estimates_sum * abundMod_results$w_aic
  pred_estimates_calc <- sum(pred_estimates_weighted)
  return(exp(pred_estimates_calc))
}

if(is_verbose){print(paste(Sys.time(), "2. finished dopar loop"))}

# is this correct -> ????
pred_per_cell <- as.data.frame(cbind(predictors_grid$id, pred_per_cell))
names(pred_per_cell) <- c("id", "abundance_pred")

# exclude NAN
pred_per_cell <- pred_per_cell[!is.nan(pred_per_cell$abundance_pred), ]

if (!is.na(exclude_year)){paste("year_to_exclude:", exclude_year)}

print(paste(Sys.time(), "sum predicted for ", year_to_predict,
            sum(pred_per_cell$abundance_pred)))
print(paste(Sys.time(),
            "range predicted for ", year_to_predict,
            range(pred_per_cell$abundance_pred)))
print(paste(Sys.time(),
            "density predicted for ", year_to_predict,
            mean(pred_per_cell$abundance_pred)))
print(paste(Sys.time(),
            "nr of values over 10 ", year_to_predict,
            sum(pred_per_cell$abundance_pred > 10)))


if (is.na(exclude_year)){
  name_suffix <- ""} else {
    name_suffix <- paste0(exclude_year, "_")
  }


saveRDS(pred_per_cell,
        file = file.path(outdir,
                         paste0("abundance_pred_per_cell_",
                                year_to_predict,"_",
                                name_suffix,
                                Sys.Date(), ".rds")))


#-----------------------#
# convert output to map #
#-----------------------#
if(is_verbose){print(paste(Sys.time(), "3. Start making map"))}


geography_for_join <- dplyr::select(geography, id, x_start, y_start)

pred_per_cell_sp <- left_join(geography_for_join,
                              pred_per_cell, by = "id") %>%
  as.data.frame()


pred_per_cell_sp <- SpatialPointsDataFrame(coords =
                                             cbind(pred_per_cell_sp$x_start,
                                                   pred_per_cell_sp$y_start),
                                           data = pred_per_cell_sp,
                                           proj4string = CRS(crs_aea),  match.ID = T)

writeOGR(pred_per_cell_sp, dsn = outdir,
         layer = paste0("prediction_shp_",
                        name_suffix,
                        year_to_predict, "_",
                        Sys.Date()), driver = "ESRI Shapefile")
# pay attention that the tif gets a new name because later we delete all
# with the same name as the shapefile
system_string <- paste0("gdal_rasterize -ot 'Float32' -a abndnc_ ",
                        "-tr 1000.0 1000.0 -l prediction_shp_",
                        name_suffix, year_to_predict, "_", Sys.Date(),
                        " ", outdir, "/",
                        "prediction_shp_", name_suffix,
                        year_to_predict, "_",
                        Sys.Date(), ".shp",
                        " ", outdir, "/",
                        "prediction_map_",name_suffix,
                        year_to_predict, "_", Sys.Date(), ".tif")
system(system_string)
#delete the shapefile again because clogs up my system
system_string_del <- paste0("rm -f ", outdir, "/prediction_shp_",
                            name_suffix, year_to_predict, "_",
                            Sys.Date(), ".*")
system(system_string_del)



save.image(file.path(outdir, paste0("abundance_pred_image_", name_suffix,
                                    year_to_predict, "_", Sys.Date(), ".RData")))

print(paste("Finish model_prediction script for year", year_to_predict,
            Sys.time()))


