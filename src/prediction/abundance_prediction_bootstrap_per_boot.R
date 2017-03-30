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
suppressPackageStartupMessages(library(stringr))  
#-----------------------------#
# command line option parsing #
#-----------------------------#
print(paste("Start abundance_prediction script", Sys.time()))


option_list <- list (
  make_option(c("-i", "--input-directory"),  dest = "input_directory",
              type = "character", help = "directory with input files",
              metavar = "/path/to/input-dir"),
  make_option("--input-directory-boot",  dest = "input_directory_boot",
                 type = "character", help = "directory with input files",
                 metavar = "/path/to/input-dir"),
  make_option(c("-o", "--output-directory"), dest = "output_directory",
               type = "character", help = "directory with output files",
               metavar = "/path/to/output-dir"),
  make_option("--nr-boot",
              dest = "nr_boot",
              type = "integer",
              help = "nr of bootstrap that is being used to predict abundance",
              metavar = "1"),
  make_option("--worker",
              dest = "worker",
              type = "integer",
              help = "number of parallel workers",
              metavar = "8"),
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

if (is.null(options$input_directory_boot)) {
  stop("input directory bootstrap not specified, check --help")
}


if (is.null(options$output_directory)) {
  stop("output directory not specified, check --help")
}


if (options$worker > 1) {
  cl <- makeCluster(options$worker, outfile = "")
  suppressMessages(registerDoParallel(cl))
}

# is quiet?
is_verbose <- options$verbose_script
# input directory
indir <- options$input_directory
if(is_verbose){print(paste("indir", indir))}

indir_boots <- options$input_directory_boot
if(is_verbose){print(paste("indir_boots", indir_boots))}


# directory in which output is written
outdir <- options$output_directory
if(is_verbose){print(paste("outdir", outdir))}

nr_boot <- as.numeric(options$nr_boot )
print(paste("nr of boots is" , nr_boot))


if(is_verbose){print(paste(Sys.time(), "0. start run"))}
#------------------------#
# command line arguments #
#------------------------#

indir_fun <- "~/orangutan_density_distribution/src/functions/"

crs_aea <- "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0
  +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # final projection


#----------------#
# Load functions #
#----------------#
# Load functions
source(file.path(indir_fun, "generic/path.to.current.R"))
source(file.path(indir_fun, "project_functions/scale.predictors.R"))


#----------------------#
# Prepare coefficients #
#----------------------#
# Load bootstrapped coefficients

all_boots_path <- path.to.current(indir_boots, "all_boots", "rds" )
if(is_verbose){print(paste("all_boots_path:", all_boots_path))}
all_boots <- readRDS(all_boots_path) %>%
  as.data.frame

#----------------------------#
# Load and prepare estimates #
#----------------------------#

# HERE it already starts to be year specific

years_to_predict <- 1999:2015

pred_per_boot <- foreach(year_to_predict = years_to_predict, .combine = rbind)  %do% {

suppressPackageStartupMessages(library(stringr))
    
predictors_grid_path <- path.to.current(indir, paste0("predictors_grid_scaled_",
                                                             year_to_predict),
                                        "rds")
if(is_verbose){print(paste("this is predictors_grid_path", predictors_grid_path))}
predictors_grid <- readRDS(predictors_grid_path)

 # predictors used in model
predictor_names <- c("year", "temp_mean", "rain_var", "rain_dry", "dom_T_OC",
                        "peatswamp", "lowland_forest",
                        "lower_montane_forest", "deforestation_hansen",
                        "human_pop_dens", "ou_killing_prediction",
                                             "perc_muslim" )

quadratic_terms_names <- c("rain_dry")

intercept <- rep(1, nrow(predictors_grid))
predictor_estimates <- cbind( intercept,
                              predictors_grid[ , predictor_names],
			                                    predictors_grid[ ,quadratic_terms_names] * predictors_grid[ ,quadratic_terms_names])

names(predictor_estimates) <- c("intercept", predictor_names,
                                paste0("I(", quadratic_terms_names, "^2)"))


#--------------------------#
# PREDICTION FOR each year #
#--------------------------#

# Alternative here to loop through id, which is added to predictor_estimates
# more correct in terms of being sure that the right thing is done,
# but takes a bit longer (52s, to 35s for 100 rows)
## PLUS PAY ATTENTION, IF PREDICTIONS NOT SAME NROW--> VALUES GET RECYCLED

if(is_verbose){print(paste("1. start pred_per_cell", Sys.time()))}

#  pred_per_cell <- foreach(cell = 1:nrow(predictor_estimates), .combine = c)  %dopar% {
  pred_per_cell <- foreach(cell = 1:700, .combine = cbind)  %dopar% {
    t_predictor_estimates <- t( predictor_estimates[cell, ])
    pred_estimates_wcoeffs  <- data.frame(mapply(`*`, all_boots[nr_boot, ], t_predictor_estimates, SIMPLIFY = F))
    pred_estimate_cell <- apply(pred_estimates_wcoeffs, 1, sum)
    pred_estimate_cell <- exp(pred_estimate_cell)
  return(pred_estimate_cell)
    }
    pred_per_boot <- c(sum(pred_per_cell),
                       min(pred_per_cell),
                       max(pred_per_cell))      
  return(pred_per_boot)
}



if(is_verbose){print(paste(Sys.time(), "2. finished dopar loop"))}

pred_per_boot <- as.data.frame(pred_per_boot)

sum_years <- as.data.frame(cbind(years_to_predict, pred_per_boot$V1))
names(sum_years) <- c("year", "sum")
min_years <- as.data.frame(cbind(years_to_predict, pred_per_boot$V2))
names(min_years) <- c("year", "min")
max_years <- as.data.frame(cbind(years_to_predict, pred_per_boot$V3))
names(max_years) <- c("year", "max")



saveRDS(sum_years,
        file = file.path(outdir,
        paste0("bootstrapped_abundance_sum_boot_",
                                nr_boot, "_",
                                Sys.Date(), ".rds")))
saveRDS(min_years,
        file = file.path(outdir,
	        paste0("bootstrapped_abundance_min_boot_",
		                                nr_boot, "_",
						Sys.Date(), ".rds")))

saveRDS(max_years,
        file = file.path(outdir,
	        paste0("bootstrapped_abundance_max_boot_",
		                                nr_boot, "_",
						Sys.Date(), ".rds")))



print(paste("Finish model_prediction script for boot", nr_boot,
            Sys.time()))


