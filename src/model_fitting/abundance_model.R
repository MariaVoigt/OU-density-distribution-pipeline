# Fit abundance model      #
#--------------------------#

#----------------#
# Load Libraries #
#----------------#

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(optparse))


#-----------------------------#
# command line option parsing #
#-----------------------------#
print(paste("Start model_fitting script", Sys.time()))


option_list <- list (
  make_option(c("-i", "--input-directory"),  dest = "input_directory",
              type = "character", help = "directory with input files",
              metavar = "/path/to/input-dir"),
  make_option(c("-o", "--output-directory"), dest = "output_directory",
              type = "character", help = "directory with output files",
              metavar = "/path/to/output-dir"),
 make_option("--ESW-aerial",
              dest = "ESW_aerial",
	      type = "double",
	      help = "aerial effective strip width"),
  make_option("--include-aerial",
              dest = "include_aerial", action="store_true",
              default=FALSE,      help="include aerial transects"),
  make_option("--stability",       action="store_true", default=FALSE,
              help="do stability analysis"),
 make_option("--exclude-year",    dest = "exclude_year", type = "integer",
             default = NA, help = "year to exclude", metavar = "2015"),
 make_option("--exclude-grid",    dest = "exclude_grid", type = "integer",
             default = NA, help = "index of grid-cells to exclude", metavar = "1"),
 make_option("--exclude-random",    dest = "exclude_rand",
             type = "integer",
             default = NA,
             help = "iteration of excluding the given percent randomly",
             metavar = "1"),
  make_option("--exclude-random-percent",
               dest = "exclude_rand_perc",
               type = "integer",
	       default = NA,
	       help = "percent cells to be excluded",
	       metavar = "10"),
  make_option(c("-q", "--quiet"), dest = "verbose_script",
              action = "store_false",
              default = TRUE,
              help = "don't print all intermediate results")
)
# verbose option a bit counterintuitive
# because I make action store_false, when I say -q that
# means that verbose == F, which is quiet

options <- parse_args(OptionParser(option_list=option_list))

if (is.null(options$input_directory)) {
  stop("input directory not specified, check --help")
}

if (is.null(options$output_directory)) {
  stop("output directory not specified, check --help")
}

exclude_year_possibilities <- c(1999:2015)

if (!is.na(options$exclude_year) && !(options$exclude_year %in% exclude_year_possibilities)) {
  stop(paste("exclude year must be between", min(exclude_year_possibilities), "and", max(exclude_year_possibilities)))
}



if (!is.na(options$exclude_rand) & is.na(options$exclude_rand_perc)){
stop(paste("exclude-grid-random needs the percent cells that have to be excluded (exclude-grid-random-percent)"))}


is_verbose <- options$verbose_script

if(is_verbose){print(paste("Reminder: Include aerial is ", options$include_aerial, ". Please
            pay attention that the same is the case for post-processing."))}
# input directory
indir <- options$input_directory
if(is_verbose){print(paste("indir", indir))}

# directory in which output is written
outdir <- options$output_directory
if(is_verbose){print(paste("outdir", outdir))}

ESW_aerial <- options$ESW_aerial
if(is_verbose){print(paste("Aerial ESW", ESW_aerial))}

do_stability <- options$stability
if(is_verbose){print(paste("stability", do_stability))}

include_aerial <- options$include_aerial
if(is_verbose){print(paste("include_aerial", include_aerial))}

exclude_year <- options$exclude_year
if(is_verbose){print(paste("exclude year", exclude_year))}

exclude_grid <- options$exclude_grid
if(is_verbose){print(paste("exclude grid", exclude_grid))}

exclude_rand <- options$exclude_rand
if(is_verbose){print(paste("exclude grid random", exclude_rand))}

exclude_rand_perc <- options$exclude_rand_perc
if(is_verbose){print(paste("exclude_rand_perc", exclude_rand_perc))}

#---------#
# Globals #
#---------#

indir_fun <- "~/orangutan_density_distribution/src/functions"
if(is_verbose){print(paste("indir_fun", indir_fun))}

cl <- makeForkCluster(outfile = "")
registerDoParallel(cl)

source(file.path(indir_fun, "project_functions/scale.predictors.R"))
source(file.path(indir_fun, "roger_functions/rogers_model_functions.R"))
source(file.path(indir_fun, "generic/path.to.current.R"))
source(file.path(indir_fun, "roger_functions/aic_c_fac.r"))
source(file.path(indir_fun, "roger_functions/get_conf_set.r"))

#define offset ground
ESW <- 0.01595  #effective strip width in km
# ESW_aerial already defined
print(paste("this is ESW aerial:", ESW_aerial))
NCS <- 1.12   #nest construction rate from Spehar et al. 2010
PNB <- 0.88  #  proportion of nest builders from Spehar et al. 2010


options("scipen" = 100, "digits" = 4)


if (is.na(exclude_year) & is.na(exclude_grid) & is.na(exclude_rand)){
  name_suffix <- ""}
   if(!is.na(exclude_year)){
     name_suffix <- paste0("year_", exclude_year, "_")}
   if(!is.na(exclude_grid)){
     name_suffix <- paste0("gridcell_", exclude_grid, "_")}
if(!is.na(exclude_rand)){
     name_suffix <- paste0("rand_", exclude_rand, "_")}

#---------------#
#  Import data  #
#---------------#

geography_path <- path.to.current(indir, "geography_observation", "rds")
if(is_verbose){print(paste("geography-path", geography_path))}
geography <- readRDS(geography_path)

transects_path <- path.to.current(indir, "transects", "rds")
if(is_verbose){print(paste("transect_path", transects_path))}
transects <- readRDS(transects_path)


predictors_path <- path.to.current(indir, "predictors_observation_20", "rds")
if(is_verbose){print(paste("predictors-path", predictors_path))}
predictors <- readRDS(predictors_path)

# calculate x and y center
geography$unscaled_x_center <-  rowMeans(cbind(geography$x_start, geography$x_end), na.rm = T)
geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = T)

#--------------------------------#
# Transform and scale predictors #
#--------------------------------#

# Transform predictors
predictors[predictors$predictor == "distance_PA", "value"] <- sqrt(
  predictors[predictors$predictor == "distance_PA", "value"])
predictors[predictors$predictor == "human_pop_dens", "value"] <- log(
  predictors[predictors$predictor == "human_pop_dens", "value"] + 1)
predictors[predictors$predictor == "deforestation_gaveau", "value"] <- sqrt(
  predictors[predictors$predictor == "deforestation_gaveau", "value"])
predictors[predictors$predictor == "plantation_distance", "value"] <- log(
  predictors[predictors$predictor == "plantation_distance", "value"] + 1)
predictors[predictors$predictor == "pulp_distance", "value"] <- log(
  predictors[predictors$predictor == "pulp_distance", "value"] + 1)
predictors[predictors$predictor == "palm_distance", "value"] <- log(
  predictors[predictors$predictor == "palm_distance", "value"] + 1)

# STARTING THE SCALING
# SCALE PREDICTORS
# these are the predictors that will be used in the model
predictor_names_for_scaling <- c( "dem",
                                 "slope",
                                 "temp_mean",
                                 "rain_dry",
                                 "rain_var",
                                 "ou_killings",
                                 "ou_killing_prediction",
                                 "human_pop_dens",
                                 "perc_muslim",
                                 "peatswamp",
                                 "lowland_forest",
                                 "lower_montane_forest" ,
                                 "road_dens",
                                 "distance_PA",
                                 "fire_dens",
                                 "deforestation_hansen",
                                 "deforestation_gaveau",
                                 "plantation_distance",
                                 "pulp_distance",
                                 "palm_distance",
                                 "plantation_age",
                                 "plantation_cover",
                                 "IOPP_age",
                                 "IOPP_cover",
                                 "ITP_age",
                                 "ITP_cover",
                                 "dom_T_OC",
                                 "dom_T_PH")

# additional predictors that have to be scaled: year and x- and y-center
predictor_names_add <- c("year", "x_center", "y_center")


# prepare predictors data-frame
predictors <- dplyr::select(predictors, id, predictor, unscaled_year = year,
                            unscaled_value = value)

# need to get rid of occurrence data
predictors <- predictors %>%
  inner_join(transects, by = "id")


# exclude aerial if that is needed
if (include_aerial == F){
  predictors <- filter(predictors, group != "aerial")
}

# delete all rows that have zero
if(is_verbose){print("how many rows with na in scaled_value")
  nrow(predictors[is.na(predictors$unscaled_value),  ])}
# deleting is.na values here
predictors <- predictors[!is.na(predictors$unscaled_value), ]

# scale predictors
predictors_obs <- scale.predictors.observation(predictor_names_for_scaling,
                                   predictor_names_add,
                                   predictors,
                                   geography)

predictors_obs <- geography %>%
  dplyr::select(-c(year,
                   unscaled_x_center,
                   unscaled_y_center)) %>%
  dplyr::select(-group) %>%
  inner_join(transects, by = "id") %>%
  inner_join(predictors_obs, by = "id")




# predictors used in model
predictor_names <- c("year", "temp_mean", "rain_var", "rain_dry", "dom_T_OC",
                     "peatswamp", "lowland_forest",
                     "lower_montane_forest", "deforestation_hansen",
                     "human_pop_dens", "ou_killing_prediction",
                     "perc_muslim" )


#  ou density and offset term for ground and absence transects
other_predictors_obs <- filter(predictors_obs, group != "aerial")
# density
other_predictors_obs$ou_dens <- (other_predictors_obs$nr_nests/
                                   (other_predictors_obs$length_km * ESW * 2))  *
  (1/(other_predictors_obs$nest_decay * NCS * PNB))
# offset
other_predictors_obs$offset_term <- log(other_predictors_obs$length_km * ESW *
                                          2 * other_predictors_obs$nest_decay *
                                          NCS * PNB)


#  ou density and offset term for aerial transects

if (include_aerial == T){
aerial_predictors_obs <- dplyr::filter(predictors_obs, group == "aerial")
# density
aerial_predictors_obs$ou_dens <- (aerial_predictors_obs$nr_nests/
                                    (aerial_predictors_obs$length_km * ESW_aerial * 2))  *
  (1/(aerial_predictors_obs$nest_decay * NCS * PNB))
# offset
aerial_predictors_obs$offset_term <- log(aerial_predictors_obs$length_km * ESW_aerial *
                                           2 * aerial_predictors_obs$nest_decay *
                                           NCS * PNB)

predictors_obs <- aerial_predictors_obs %>%
  bind_rows(other_predictors_obs)
}else{
  predictors_obs <- other_predictors_obs
}

# bind the two together
predictors_obs <- predictors_obs %>%
  arrange(id) %>%
  as.data.frame(.)


if(is_verbose){print("look at predictors_obs")
  str(predictors_obs)
  summary(predictors_obs)}

# save the relevant output for the prediction and the validation (USED IN THESE SCRIPTS)
saveRDS(predictors_obs, file = file.path(outdir, paste0("predictors_observation_scaled_",
                                                          name_suffix,
                                                          Sys.Date(), ".rds")))


# now exclude the year that needs to be excluded
if (!is.na(exclude_year)){
    predictors_excluded <- predictors_obs[predictors_obs$unscaled_year == exclude_year, ]
    predictors_obs <- predictors_obs[predictors_obs$unscaled_year != exclude_year, ]
    nr_excluded <- nrow(predictors_excluded)}

# or the grid_cell
if (!is.na(exclude_grid)){
  predictors_exclude <- predictors_obs[predictors_obs$grid_id == exclude_grid, ]
  predictors_obs <- predictors_obs[predictors_obs$grid_id != exclude_grid, ]
  nr_excluded <- nrow(predictors_excluded)
}

if (!is.na(exclude_rand)){
  ids_to_exclude <- sample(predictors_obs$id,
                           size = nrow(predictors_obs)/100 * exclude_rand_perc,
                           replace = FALSE)
  predictors_excluded <- predictors_obs[predictors_obs$id %in% ids_to_exclude, ]
  predictors_obs <- predictors_obs[!predictors_obs$id %in% ids_to_exclude, ]
  nr_excluded <- nrow(predictors_excluded)
}

# also we increase maxit for the two cases,
# because then slightly less data
if (is.na(exclude_grid) & is.na(exclude_year) & is.na(exclude_rand)){
    nr_maxit <- 250}else{
                   nr_maxit <- 500}
if(is_verbose){ print(paste("3. start making all_model_terms", Sys.time()))}

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
                                         "deforestation_hansen",
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
             "deforestation_hansen",
             "human_pop_dens",
             "ou_killing_prediction",
             "perc_muslim",
             "I(rain_dry^2)")


# save model_terms here
model_terms <- names(glm.nb(as.formula(paste("nr_nests~", paste(m_terms,
                                                                collapse = "+"),
                                             "+ offset(offset_term)",
                                             sep = "")),
                              data = predictors_obs,
                            control = glm.control(maxit = nr_maxit))$coefficients)


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
  if(is_verbose){print(paste("Start stability calculation", Sys.time()))}
full_model <- paste(
  m_terms[all_model_terms[nrow(all_model_terms), ] == 1],
  collapse = "+")
if(is_verbose){print(paste("This is the full-model", full_model))}
model <- as.formula(
  paste("nr_nests ~", full_model, "+ offset(offset_term)"))

res_full <- glm.nb(model, data = predictors_obs,

                   control = glm.control(maxit = nr_maxit))


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
if(is_verbose){print(paste("8. Start running models", Sys.time()))}

results_res <- foreach(i = 1:nrow(all_model_terms),
                       .combine = rbind) %dopar%{
    # make results dataframe
    if (is.na(exclude_year) & is.na(exclude_grid) & is.na(exclude_rand)){
      result <- as.data.frame(matrix(NA, ncol = 3 *
                                       length(model_terms) + 6,
                                     nrow = 1))
       names(result) <- c("model", paste("coeff", model_terms, sep = "_"),
                                              paste("P",model_terms,sep = "_"),
                                              paste("SE", model_terms, sep = "_"),
                                              "theta", "SE.theta", "AIC", "R2", "nr_excluded"
                           )} else {
      result <- as.data.frame(matrix(NA, ncol = 3 * length(model_terms) + 6,
                                                          nrow = 1))
      names(result) <- c("model", paste("coeff", model_terms, sep = "_"),
                                              paste("P",model_terms,sep = "_"),
                                              paste("SE", model_terms, sep = "_"),
                                              "theta", "SE.theta", "AIC", "R2",
                                              "R2_cross")
                         }
    # make model
    model <- as.formula(
            paste("nr_nests ~",
	    paste(m_terms[all_model_terms[i, ] == 1], collapse = "+"),
            "+ offset(offset_term)"))
    res <- glm.nb(model, data = predictors_obs,
                          control = glm.control(maxit = nr_maxit))


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
    #all_model_terms[i, ]==1
    predictors_obs_pred <- predictors_obs
    predictors_obs_pred$offset_term <- 0
    # Comparison of observed data vs prediction
    prediction_per_transect <-  predict.glm(res,
                                            newdata = predictors_obs_pred,
                                            type = "response")

     comparison_lm = lm(log(predictors_obs$ou_dens + 1) ~
                        log(prediction_per_transect + 1) )

        result[ , "R2"] <- summary(comparison_lm)$r.squared
    # if we are excluding years, this is the test of predicted data vs observed data
    # for this year (with which the model wasn't fitted)
        # probably I could code this bit similarly for all


        if (!is.na(exclude_year) | !is.na(exclude_grid) | !is.na(exclude_rand)){
           predictors_excluded_pred <- predictors_excluded
          predictors_excluded_pred$offset_term <- 0
          prediction_transect_excluded <-  predict.glm(res,
                                                            newdata = predictors_excluded_pred,
                                                            type = "response")
          cross_lm = lm(log(predictors_excluded$ou_dens+ 1) ~
                               log(prediction_transect_excluded + 1))

          result[ , "R2_cross"] <- summary(cross_lm)$r.squared
         result[ , "nr_excluded"] <- nr_excluded
         }

  return(result)
}


c_set <- cbind(as.character(results_res$model), conf.set(aic = results_res$AIC) )
names(c_set) <- c("model", names(c_set)[-1])

# for the export
abundMod_results <- results_res %>%
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
saveRDS(abundMod_results, file = file.path(outdir, paste0("abundMod_results_",
                                                          name_suffix,
                                                Sys.Date(), ".rds")))

# these are the terms that go into the model (need to guarantee same things in validation)
saveRDS(m_terms, file = file.path(outdir, paste0("m_terms_",
                                                 name_suffix,
                                                 Sys.Date(), ".rds")))

# save the model results for interpretation
write.csv(results_out,
          file = file.path(outdir,
                           paste0("abundMod_results_",
                                  name_suffix,
                                  Sys.Date(), ".csv")))
# save the mean coefficients for interpretation
write.csv(summary_mean_coefficients,
          file = file.path(outdir,
                           paste0("abundMod_mean_coefficients_",
                                  name_suffix,
                                  Sys.Date(), ".csv")))


#save.image(file.path(outdir, paste0("abundance_model_fitting_",
#                                    name_suffix,
# Sys.Date(), ".RData")))

print(paste("Finished model_fitting script, at", Sys.time()))
