#--------------------------#
# Fit abundance model      #
#--------------------------#
rm(list = ls())
gc()


#----------------#
# Load Libraries #
#----------------#
library(dplyr)
library(tidyr)
library(stringr)
library(MASS)
library(gtools)
library(reshape2)
library(parallel)
library(foreach)
library(doParallel)
library(pscl)


#---------#
# Globals #
#---------#

args <- commandArgs(trailingOnly = TRUE)
print(paste("args", args))
# the directory in which the files are

indir <- args[1]

print(paste("indir ", indir))
# directory in which output is written

outdir <- args[2]
print(paste("outdir ", outdir))

do_stability <- args[3]

if(is.na(do_stability)){do_stability <- FALSE} else{
  if (do_stability == "do_stability"){do_stability <- TRUE} else{
    stop("there is the wrong input in do_stability")
  }
}

indir_fun <- "../functions"
print(paste("indir_fun", indir_fun))

cl <- makeForkCluster(outfile = "")
registerDoParallel(cl)

source(file.path(indir_fun, "roger_functions/rogers_model_functions.R"))
source(file.path(indir_fun, "generic/path.to.current.R"))
source(file.path(indir_fun, "roger_functions/aic_c_fac.r"))
source(file.path(indir_fun, "roger_functions/get_conf_set.r"))

#define offset ground
ESW <- 0.01595  #effective strip width in km
NCS <- 1.12   #nest construction rate from Spehar et al. 2010
PNB <- 0.88  #  proportion of nest builders from Spehar et al. 2010


options("scipen" = 100, "digits" = 4)


#---------------#
#  Import data  #
#---------------#

geography_path <- path.to.current(indir, "geography_observation", "rds")
print(paste("geography-path", geography_path))
geography <- readRDS(geography_path)

transects_path <- path.to.current(indir, "transects", "rds")
print(paste("transect_path", transects_path))
transects <- readRDS(transects_path)


predictors_path <- path.to.current(indir, "predictors_observation", "rds")
print(paste("predictors-path", predictors_path))
predictors <- readRDS(predictors_path)


# these are the predictors that will be used in the model
predictor_names <- c("year", "temp_mean", "rain_var", "rain_dry", "dom_T_OC",
                     "peatswamp", "lowland_forest",
                     "lower_montane_forest", "deforestation",
                     "human_pop_dens", "ou_killing_prediction",
                     "perc_muslim" )

geography <- dplyr::select(geography, -year)

print("how many rows with na in scaled_value")
nrow(predictors[is.na(predictors$scaled_value),  ])
# deleting is.na values here
predictors <- predictors[!is.na(predictors$scaled_value), ]
predictors_obs <- predictors %>%
  dplyr::filter(predictor %in% predictor_names) %>%
  dcast(id + year ~ predictor,  value.var = "scaled_value")%>%
  inner_join(geography, by = "id")%>%
  dplyr::select(-group) %>%
  inner_join(transects, by = "id" ) %>%
  dplyr::filter(group != "aerial")


predictors_obs$ou_dens <- (predictors_obs$nr_nests/ (predictors_obs$length_km * ESW * 2))  *
  (1/(predictors_obs$nest_decay * NCS * PNB))
# HERE WE ARE ONLY LOOKING AT NESTS
predictors_obs$offset_term <- log(predictors_obs$length_km * ESW * 2 * predictors_obs$nest_decay * NCS * PNB)


# SCALE YEAR
scaled_year <- as.vector(scale(predictors_obs$year))
predictors_obs$unscaled_year <- predictors_obs$year
predictors_obs$year <- as.numeric(scaled_year)

# calculate x and y center
predictors_obs$x_center <-  rowMeans(cbind(predictors_obs$x_start, predictors_obs$x_end), na.rm = T)
predictors_obs$y_center <- rowMeans(cbind(predictors_obs$y_start, predictors_obs$y_end), na.rm = T)

# calculate ou_dens
predictors_obs$nr_ou_per_km2 <- predictors_obs$nr_nests /
  (predictors_obs$length_km * ESW * 2 * predictors_obs$nest_decay  * NCS * PNB )

print("look at predictors_obs")
str(predictors_obs)
summary(predictors_obs)

 print(paste("3. start making all_model_terms", Sys.time()))

 # #build models needed for analysis with a function
all_model_terms <- built.all.models(env.cov.names =
                                      c( "year",
                                         "temp_mean",
                                         "rain_var",
                                         "rain_dry",
                                         "dom_T_OC",
                                         "peatswamp",
                                         "lowland_forest",
                                         "lower_montane_forest",
                                         "deforestation",
                                         "human_pop_dens",
                                         "ou_killing_prediction",
                                         "perc_muslim"),
                                    env.cov.int = list(),
                                    env.cov.2 = c("rain_dry"))



 print(paste("4. end making all model terms", Sys.time()))


m_terms <- c("1",
             "year",
             "temp_mean",
             "rain_var",
             "rain_dry",
             "dom_T_OC",
             "peatswamp",
             "lowland_forest",
             "lower_montane_forest",
             "deforestation",
             "human_pop_dens",
             "ou_killing_prediction",
             "perc_muslim",
             "I(rain_dry^2)")


# save model_terms here
model_terms <- names(glm.nb(as.formula(paste("nr_nests~", paste(m_terms,
                                                                collapse = "+"),
                                             "+ offset(offset_term)",
                                             sep = "")),
                              data = predictors_obs)$coefficients)
# prediction estimates
intercept <- rep(1, nrow(predictors_obs))
predictor_estimates <- cbind( intercept,
                               predictors_obs[ , predictor_names],
                               predictors_obs[ ,"rain_dry"] *
                                 predictors_obs[ ,"rain_dry"])

names(predictor_estimates) <- c("intercept", predictor_names,
                                paste0("I(", "rain_dry", "^2)"))



# calculate stability of the full model if desired
if(do_stability){
print(paste("Start stability calculation", Sys.time()))
full_model <- paste(
  m_terms[all_model_terms[nrow(all_model_terms), ] == 1],
  collapse = "+")
print(paste("This is the full-model", full_model))
model <- as.formula(
  paste("nr_nests ~", full_model, "+ offset(offset_term)"))

res_full <- glm.nb(model, data = predictors_obs)

# HERE I CAN NOW USE THE OTHER FUNCTION
dfbeta_frame <- data.frame(slope=res_full$coefficients, res_full$coefficients+
                             t(apply(X=dfbeta(res_full),
                                     MARGIN=2, FUN=range)))
names(dfbeta_frame) <- c("slope", "min", "max")

write.csv(dfbeta_frame, file.path(outdir,
                                  paste0("glm_abundance_stability_",
                                         Sys.Date(), ".csv")), row.names = T)
}

# #run models
print(paste("8. Start running models", Sys.time()))

results_res <- foreach(i = 1:nrow(all_model_terms), .combine = rbind) %dopar% {
#system.time(results_res <- foreach (i = 1:10, .combine = rbind) %do% {
    # create objects for storing results
    # modelinfo + schätzung von coeffizienten und pi values
    result <- as.data.frame(matrix(NA, ncol = 3 * length(model_terms) + 5,
                                        nrow = 1))
    names(result) <- c("model", paste("coeff", model_terms, sep = "_"),
                            paste("P",model_terms,sep = "_"),
                            paste("SE", model_terms, sep = "_"),
                       "theta", "SE.theta", "AIC", "R2")

    # model fitting
    model <- as.formula(
        paste("nr_nests ~",
              paste(m_terms[all_model_terms[i, ] == 1], collapse = "+"),
              "+ offset(offset_term)"))
    res <- glm.nb(model, data = predictors_obs)

    # model
    result[ , "model"] <- paste(m_terms[all_model_terms[i, ] == 1], collapse = "+")


    # coefficients
    model_coefficients <- as.vector(res$coefficients)
    result[ , paste0("coeff_", names(res$coefficients))] <- model_coefficients


    # p value
    result[ , paste0("P_", names(res$coefficients))] <-
      summary(res)$coefficients[ , "Pr(>|z|)"][names(
        res$coefficients)] #w/o parameter for autocor.

    # SE
    result[ , paste0("SE_", names(res$coefficients))] <-
      summary(res)$coefficients[ , "Std. Error"][names(
        res$coefficients)]#add line for SE

    # theta
   result[ , "theta"] <- summary(res)$theta
   result[ , "SE.theta"] <- summary(res)$SE.theta

    # aic in last column,
    result[ , "AIC"] <- extractAIC(res)[2]

    # what do I need to do
    # I need to get only the prediction estimates columns that are true in
    # all_model_terms[i, ]==1
    predictors_obs_pred <- predictors_obs
    predictors_obs_pred$offset_term <- 0
    # prediction estimates
    prediction_per_transect <-  predict.glm(res,
                                            newdata = predictors_obs_pred,
                                            type = "response")

    comparison_lm = lm(predictors_obs$nr_ou_per_km2 ~ prediction_per_transect )

    result[ , "R2"] <- summary(comparison_lm)$r.squared
    return(result)
}


c_set <- cbind(as.character(results_res$model), conf.set(aic = results_res$AIC) )
names(c_set) <- c("model", names(c_set)[-1])

# for the export
abundMod_result <- results_res %>%
  mutate(w_aic = c_set$w.aic)


# calculate coefficient summary
summary_mean_coefficients <- calculate.mean.coefficients(m_terms, results_res, c_set)

# make table for model output
c_set$model <- NULL
names(c_set)[1] <- "AIC"
names(c_set)
c_set$d.aic <- NULL
c_set$w.aic <- NULL
results_out <- right_join(results_res, c_set,  by="AIC")
results_out <- results_out[order(results_out$AIC), ]


# save the relevant output for the prediction and the validation
saveRDS(abundMod_result, file = file.path(outdir, paste0("abundMod_results",
                                                Sys.Date(), ".rds")))

# these are the terms that go into the model (need to guarantee same things in validation)
saveRDS(m_terms, file = file.path(outdir, paste0("m_terms",
                                                 Sys.Date(), ".rds")))

# save the model results for interpretation
write.csv(results_out,
          file = file.path(outdir,
                           paste0("abundMod_results",
                                               Sys.Date(), ".csv")))
# save the mean coefficients for interpretation
write.csv(summary_mean_coefficients,
          file = file.path(outdir,
                           paste0("abundMod_mean_coefficients",
                                  Sys.Date(), ".csv")))


save.image(file.path(outdir, paste0("abundance_model_fitting_", Sys.Date(), ".RData")))
print(paste("11. finished script, finally, at", Sys.time()))
