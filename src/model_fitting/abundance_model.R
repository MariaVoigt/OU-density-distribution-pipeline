#---------------------#
# Fit abundance model #
#---------------------#
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
# to check whether it is right (in log-script)
print(paste("indir ", indir))
# directory in which output is written
outdir <- args[2]
print(paste("outdir ", outdir))

indir_fun <- "../functions"
print(paste("indir_fun", indir_fun))


cl <- makeForkCluster(outfile = "")
registerDoParallel(cl)

source(file.path(indir_fun, "/roger_functions/rogers_model_functions.R"))
source(file.path(indir_fun, "/generic/path.to.current.R"))
source(file.path(indir_fun, "/roger_functions/aic_c_fac.r"))
source(file.path(indir_fun, "/roger_functions/get_conf_set.r"))

#define offset ground
ESW <- 0.01595  #effective strip width in km
NCS <- 1.12   #nest construction rate from Spehar et al. 2010
PNB <- 0.88  #  proportion of nest builders from Spehar et al. 2010


options("scipen" = 100, "digits" = 4)


#---------------#
#  Import data  #
#---------------#

geography_path <- path.to.current(indir, "geography_observation", "rds")
geography <- readRDS(geography_path)

transects_path <- path.to.current(indir, "transects", "rds")
transects <- readRDS(transects_path)

predictors_path <- path.to.current(indir, "predictors_observation", "rds")
predictors <- readRDS(predictors_path)


# these are the predictors that will be used in the model
predictor_names <- c("year", "temp_mean", "rain_var", "rain_dry", "dom_T_OC",
                     "dom_T_PH", "peatswamp", "lowland_forest",
                     "lower_montane_forest", "deforestation", "fire_dens",
                     "road_dens", "distance_PA", "human_pop_dens",
                     "ou_killing_prediction", "perc_muslim" )

geography <- dplyr::select(geography, -year)

data <- predictors %>%
  dplyr::filter(predictor %in% predictor_names) %>%
  dcast(id + year ~ predictor,  value.var = "scaled_value")%>%
  inner_join(geography, by = "id")%>%
  dplyr::select(-group) %>%
  inner_join(transects, by = "id" )

# work on aerial transects
# here we need to calculate first the aerial index (nests / km) and then transform it to nest-density
# then we include the offset without length_km * 2 * ESW
# for the transect this goes into the term
aerial_data <- dplyr::filter(data, group == "aerial")
# Ai is aerial index (number of nests detected per kilometer of flight)
aerial_data$AI <- aerial_data$nr_nests / aerial_data$length_km
# calculate orangutan nest density from aerial index with formula 6 given in Ancrenaz et al., 2004
aerial_data$nr_nests <- round(exp(4.7297 + 0.9796 * log(aerial_data$AI)))

aerial_data$ou_dens <- aerial_data$nr_nests / (aerial_data$nest_decay * NCS * PNB)
aerial_data$offset_term <- log(1 * 1* aerial_data$nest_decay * NCS * PNB )

other_data <- filter(data, group != "aerial")
other_data$ou_dens <- (other_data$nr_nests/ (other_data$length_km * ESW * 2))  *
  (1/(other_data$nest_decay * NCS * PNB))

other_data$offset_term <- log(other_data$length_km * ESW * 2 * other_data$nest_decay * NCS * PNB)
names_data <- names(other_data)

data <- aerial_data %>%
  dplyr::select(id:length_km, nr_nests, nest_decay, ou_dens, offset_term)
print("This has to be true:")
unique(names(data) == names(other_data))
# HAS TO BE TRUE
data <- data %>%
  bind_rows(other_data) %>%
  arrange(id)


data$z.year <- scale(data$year)
data$unscaled_year <- as.numeric(data$year)
data$year <- as.numeric(data$z.year)
data$z.year <- NULL


# calculate x and y center
data$x_center <-  rowMeans(cbind(data$x_start, data$x_end), na.rm = T)
data$y_center <- rowMeans(cbind(data$y_start, data$y_end), na.rm = T)

print("look at data")
str(data)
summary(data)

 print(paste("3. start making all_model_terms", Sys.time()))

 # #build models needed for analysis with a function
all_model_terms <- built.all.models(env.cov.names =
                                      c( "year",
                                         "temp_mean",
                                         "rain_var",
                                         "rain_dry",
                                         "dom_T_OC",
                                         "dom_T_PH",
                                         "peatswamp",
                                         "lowland_forest",
                                         "lower_montane_forest",
                                         "deforestation",
                                         "fire_dens",
                                         "road_dens",
                                         "distance_PA",
                                         "human_pop_dens",
                                         "ou_killing_prediction",
                                         "perc_muslim"),
                                    env.cov.int = list(c("year", "deforestation"),
                                                       c("year", "road_dens"),
                                                       c("year", "human_pop_dens"),
                                                       c("year", "distance_PA")),
                                    env.cov.2 = c("temp_mean",
                                                  "rain_dry"))



 print(paste("4. end making all model terms", Sys.time()))


m_terms <- c("1",
             "year",
             "temp_mean",
             "rain_var",
             "rain_dry",
             "dom_T_OC",
             "dom_T_PH",
             "peatswamp",
             "lowland_forest",
             "lower_montane_forest",
             "deforestation",
             "fire_dens",
             "road_dens",
             "distance_PA",
             "human_pop_dens",
             "ou_killing_prediction",
             "perc_muslim",
             "deforestation:year",
             "road_dens:year",
             "human_pop_dens:year",
             "distance_PA:year",
             "I(temp_mean^2)",
             "I(rain_dry^2)")


# save model_terms here
model_terms <- names(zeroinfl(as.formula(paste("nr_nests~", paste(m_terms, collapse = "+"),
                                             "+ offset(offset_term) | 1", sep = "")),
                              data = data, dist = "negbin")$coefficients$count)
# offset does not appear in model_terms
model_terms <- c(model_terms ,"ac_term")


#get autocorrelation term
#run full model and derive autocorrelation term
res <- zeroinfl(as.formula(paste("nr_nests~", paste(m_terms, collapse = "+"),
                              "+ offset(offset_term) | 1", sep = "")) ,
             data = data, dist = "negbin")
resis <- residuals(res)


get_wsd <- function(xsd){
  xac_term <- get.1d.ac(resis = resis, ac.sd=xsd, lat = data$y_center,
                        long = data$x_center,
                        contr.fac=NULL)
  xac_term <- as.vector(scale(xac_term))
  xres <- zeroinfl(as.formula(paste(paste("nr_nests~", paste(m_terms, collapse = "+"),
                          "xac_term", "offset(offset_term) ", sep = "+"), "| 1",
                          sep = "")),
         data = data, dist = "negbin")
  aic <- -2 * xres$loglik + 2 *
    length(coefficients(xres)) +
    aic.c.fac(N = nrow(data),
              k = length(coefficients(xres)))
  return(aic)}

 print(paste("5. start estimating autocorrelation term", Sys.time()))

all_sd <- seq(100, 40000, by = 100)
system.time(
    all_aic <- unlist(lapply(all_sd, get_wsd))
)
warnings()

save.image(file.path(outdir, "image_post_aic.RData"))

png(file.path(outdir, 'aic_curve.png'), type = 'cairo')
plot(all_sd, all_aic)
dev.off()
print("which is the minimum aic and where is it")

all_aic[which.min(all_aic)]
all_sd[which.min(all_aic)]

 print(paste("6. exported aic optimum plot", Sys.time()))
# # where is the AIC minimum in the plot--> check the minimum in this range
w.sd <- optimize(get_wsd, lower = 10000, upper = 20000)
ac_term <- get.1d.ac(resis = resis, ac.sd = w.sd$minimum, lat = data$y_center,
                     long = data$x_center, contr.fac=NULL)
ac_term <- as.vector(scale(ac_term))

print("this is ac_term:")
str(ac_term)
print(paste("7. finished aic term", Sys.time()))
saveRDS(ac_term, file.path(outdir, paste0("ac_term_", Sys.Date(), ".rds")))
save.image(file.path(outdir, "abundance_model_ac_term_image.RData"))

# #run models
print(paste("8. start running models", Sys.time()))

results_res <- foreach(i = 1:nrow(all_model_terms), .combine = rbind,
                       .export = c("ac_term")) %dopar% {

    # create objects for storing results
    # modelinfo + schätzung von coeffizienten und pi values
    result <- as.data.frame(matrix(NA, ncol = 3 * length(model_terms) + 2,
                                        nrow = 1))
    names(result) <- c("model", paste("coeff", model_terms, sep = "_"),
                            paste("P",model_terms,sep = "_"),
                            paste("SE", model_terms, sep = "_"), "AIC")

    # check here with space before + ac_term
    model <- as.formula(
        paste("nr_nests ~",
              paste(m_terms[all_model_terms[i, ] == 1], collapse = "+"),
              "+ ac_term + offset(offset_term) | 1"))
    res <- zeroinfl(model, data = data, dist = "negbin")

    # model
    result[ , "model"] <- paste(m_terms[all_model_terms[i, ] == 1], collapse = "+")

    # coefficients (+ 1 because first column are the models)
    result[ , paste0("coeff_", names(res$coefficients$count))] <-
      as.vector(res$coefficients$count)

    # p value
    result[ , paste0("P_", names(res$coefficients$count))] <-
      summary(res)$coefficients$count[ , "Pr(>|z|)"][names(
        res$coefficients$count)] #w/o parameter for autocor.

    # SE
    result[ , paste0("SE_", names(res$coefficients$count))] <-
      summary(res)$coefficients$count[ , "Std. Error"][names(
        res$coefficients$count)]#add line for SE

    # aic in last column, roger script 10b_aic_mmi_applied_handout
    aic <- -2 * res$loglik + 2 *
      length(coefficients(res)) +
      aic.c.fac(N = nrow(data),
                k = length(coefficients(res)))

    result[ , "AIC"] <- aic

    return(result)
}

dim(results_res)
save.image(file.path(outdir, "abundance_model_ac_term_image_post_foreach_doparloop.RData"))

c_set <- cbind(as.character(results_res$model), conf.set(aic = results_res$AIC) )
names(c_set) <- c("model", names(c_set)[-1])

# for the export
export <- results_res %>%
  dplyr::select(model,
                contains("coeff")) %>%
  mutate(w_aic = c_set$w.aic)

# results_res
saveRDS(export, file = file.path(outdir, paste0("coeff_weights_", Sys.Date(), ".rds")))


#calculate summary stats
#parameter estimates
# !!! ATTENTION
# here all parameter coeff_, excluding the model, including the ac_term
length_par_est <- length(m_terms) + 1
par.est <- results_res[ , 2:(length_par_est + 1)] # +1 because we start at 2
#not considering models without parameter estimate
par.est.av.1 <- c()
#
# # here we loop through par.est, so it has to have same length
for(i in 1:length(par.est)) {
  temp.par.est <- cbind(par.est[ ,i], c_set$w.aic)
  temp.par.est <- subset(temp.par.est, !is.na(temp.par.est[ , 1]))
  temp.par.est[ , 2] <- temp.par.est[, 2] / sum(temp.par.est[ , 2])
  par.est.av.1 <- c(par.est.av.1,
                    weighted.mean(temp.par.est[ ,1],
                                  temp.par.est[ ,2],
                                  na.rm = T))
}


#including models without parameter estimate
par.est.no.NA <- results_res[ , 2:(length_par_est + 1)]

for (i in 1:length(par.est.no.NA)) {
  par.est.no.NA[,i][is.na(par.est.no.NA[,i])] <- 0
}
par.est.av.2 <- apply(par.est.no.NA * c_set$w.aic, 2, sum)
res.par.est.av <- data.frame(par.est.av.1, par.est.av.2)

print(paste("10. saving output", Sys.time()))

# # join c_set and results_res for showing
c_set$model <- NULL
names(c_set)[1] <- "AIC"
names(c_set)
results_out <- right_join(results_res, c_set,  by="AIC")
str(results_out)
results_out <- results_out[order(results_out$AIC), ]

write.csv(results_out, file = file.path(outdir, paste0("model_results_abundance_model_", Sys.Date(), ".csv")))
write.csv(res.par.est.av, file = file.path(outdir, paste0("av_parameter_estimates_abundance_model_", Sys.Date(), ".csv")))

save.image(file.path(outdir, paste0("abundance_model_and_prediction_", Sys.Date(), ".RData")))

print(paste("11. finished script, finally, at", Sys.time()))
