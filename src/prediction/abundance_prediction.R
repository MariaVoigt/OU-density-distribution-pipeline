
#-----------------------------------------------#
# Model prediction:                             #
# the next part is for making model predictions #
# throughout study area                         #
# define predictors of 'grid'                   #
#-----------------------------------------------#
#----------------#
# Load Libraries #
#----------------#

library(parallel)
library(foreach)
library(doParallel)
library(reshape2)
library(dplyr)


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
indir_fun <- args[3]
print(paste("indir fun ", indir_fun))
indir_predictors <- args[4]
print(paste("indir predictors", indir_predictors))
year_to_predict <- as.numeric(args[5])
print(paste("year " , year_to_predict))



#-----------#
# LOAD DATA #
#-----------#
# Load functions
source(file.path(indir_fun, "path.to.current.R"))
print("function loaded")

# Load coefficients and weights
coeff_weights_path <- path.to.current(indir, "coeff_weights_abundance", "rds" )
print(paste0("this is coeff_weights_path", coeff_weights_path))
coeff_weights <- readRDS(coeff_weights_path)
# exclude the first column, which contains models, exclude the coefficient of the
# autocorellation term and the weighted aic of the model
# we are excluding the autocorellation term, because the mean is zero, and thus
# it will not have an influence on the overall outcome
coeffs <- coeff_weights[ , 2:(length(coeff_weights) - 2)]
# all predictors that are not included in a specific model have NA as their
# coefficient, which is replaced by 0, so that the estimate is also 0
coeffs[is.na(coeffs) == T] <- 0

# Load estimates for observation and grid
# these are the predictors that will be used in the prediction
# THEY MUST BE IN THE ORDER IN WHICH THEY APPEAR IN THE COEFFICIENTS
# CODE THIS
# here only separate after first _
predictor_names_coeffs <- gsub("coeff_","", names(coeffs))
interaction_terms_names <- predictor_names_coeffs[predictor_names_coeffs %in%
                           paste0("year:", predictor_names_coeffs)] 
interaction_terms_names <- gsub("year:", "", interaction_terms_names)
quadratic_terms_names <- predictor_names_coeffs[predictor_names_coeffs %in%
                                                 paste0("I(", predictor_names_coeffs, "^2)")]
quadratic_terms_names <- gsub("I\\(|\\^2\\)", "", quadratic_terms_names )


predictor_names_coeffs <- predictor_names_coeffs[predictor_names_coeffs != "(Intercept)"]
# don't include interaction or quadratic term
predictor_names <- predictor_names_coeffs[!grepl("I(*)", predictor_names_coeffs)]
predictor_names <- predictor_names[!grepl("year[:punct:]*", predictor_names)]
# this is not a good fix, the problem is that with the second grepl also variable "year" gone
predictor_names <- c("year", predictor_names)
                                        # predictors for year on grid

predictors_path <- path.to.current(indir_predictors, paste0("predictors_occ_grid_",
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
                    predictors[ ,quadratic_terms_names] * predictors[ ,quadratic_terms_names],
                    predictors[ , "year"] * predictors[ , interaction_terms_names])

names(predictor_estimates) <- c("intercept", predictor_names,
                                paste0("I(", quadratic_terms_names, "^2)"),
				paste0("year:", interaction_terms_names))



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
    pred_estimates_calc <- exp(pred_estimates_sum) * 1/(1 + exp(-(pred_estimates_sum)))
                 pred_estimates_weighted <- pred_estimates_calc * coeff_weights$w_aic
                 sum(pred_estimates_weighted)
                }

print(paste(Sys.time(), "2. finished dopar loop"))


pred_per_cell <- as.data.frame(cbind(predictors$id, pred_per_cell))
names(pred_per_cell) <- c("id", "abundance_pred")

# exclude NAN
pred_per_cell <- pred_per_cell[!is.nan(pred_per_cell$abundance_pred), ]


print(paste(Sys.time(), "sum predicted for ", year_to_predict,
            sum(pred_per_cell$abundance_pred)))
print(paste(Sys.time(), "range predicted for ", year_to_predict,
            range(pred_per_cell$abundance_pred)))
# has to be between 0 and 1

save.image(file.path(outdir, paste0("abundance_pred_image_", year_to_predict, "_",
                                    Sys.Date(), ".RData")))
saveRDS(pred_per_cell,
          file = file.path(outdir,
                           paste0("abundance_pred_per_cell_",
                                  year_to_predict,"_",
                                  Sys.Date(), ".rds")))


print(paste(Sys.time(), "3. wrote results and done :-)"))
