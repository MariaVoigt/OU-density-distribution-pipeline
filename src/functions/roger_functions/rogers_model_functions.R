#######functions needed to run script
#function to build all models from given predictors and response
built.all.models <- function(env.cov.names, env.cov.int, env.cov.2){
  all.mods <- permutations(n = 2,r = length(env.cov.names) + length(env.cov.int) +
                             length(env.cov.2), v = c(0, 1), repeats.allowed = T)
  colnames(all.mods) <- c(env.cov.names,
                          unlist(lapply(env.cov.int, function(x){paste(x,collapse=":")})),
                          unlist(lapply(env.cov.2, function(x){paste(c("I(", x, "^2)"),
                                                                     collapse="")})))
  # the problem is, that we cannot include interactions, if covariates are not
  # in the model, so we can only use interaction a*b, if a and b are in the model
  all.problems <- c(env.cov.int, lapply(env.cov.2, function(x){c(x,x)}))
  if (length(all.problems) > 0){
    for(i in 1:length(all.problems)){
      to.del <- all.mods[ , length(env.cov.names) + i] == 1 &
        (all.mods[ ,colnames(all.mods) == (all.problems[[i]])[1]] == 0 |
           all.mods[ ,colnames(all.mods) == (all.problems[[i]])[2]] == 0
        )
      all.mods <- all.mods[to.del == F,     ]
    }
  }
  all.mods <- cbind("(Intercept)" = rep(1, nrow(all.mods)), all.mods)
  return(all.mods)
}

#function to calculate autocorrelation term
get.ac <- function(measure, lat, lon=NULL, c.fac=NULL, same.pos.allowed = T,
                   exclude.ident.ID = F, rowID = NULL) {
  if (length(lon) == 0) {lon = rep(0, length(measure))}
  if(length(c.fac) == 0) {c.fac = rep(1, length(measure))}#create constant control factor in case it doesn't exist
  if (length(rowID) == 0) {rowID = 1:length(measure)}
  if (same.pos.allowed == T) {same.pos.allowed=1} else{same.pos.allowed = 0}
  ac.tot = rep(NA,length(measure))
  xlevels = names(table(c.fac)) #determine levels of control factor
  for (k in 1:length(xlevels)) { #for each level of the control factor
    sel.rows <- (1:length(measure))[c.fac == xlevels[k]]
    sel.measure <- measure[sel.rows]
    sel.lat <- lat[sel.rows]
    sel.lon <- lon[sel.rows]
    selID <- rowID[sel.rows]
    ac.ind <- c()
    for (i in 1:length(sel.measure)) {
      act.lat <- sel.lat[i]; act.lon <- sel.lon[i]; actID <- selID[i]
      if (exclude.ident.ID == F) {
        ind.weight <- 1 / (same.pos.allowed + sqrt((act.lat - sel.lat)^2 +
                                                     (act.lon - sel.lon)^2))
        ac.ind <- c(ac.ind, sum(sel.measure[-i] * ind.weight[-i]) /
                      sum(ind.weight[-i]))
      }else{
        ind.weight <- same.pos.allowed + sqrt((act.lat - sel.lat)^2 +
                                                (act.lon - sel.lon)^2)
        ac.ind <- c(ac.ind, sum(sel.measure[selID != actID] /
                                  ind.weight[selID != actID]) /
                      sum(1 / ind.weight[selID != actID]))
      }
    }
    ac.tot[sel.rows] <- ac.ind
  }
  return(ac.tot)
}

###function to interpolate autocorrelation term across study area
predict.ac <- function(ac_term, lat, lon = NULL, same.pos.allowed = T,
                       predict.lat, predict.lon){
  if (length(lon) == 0){lon <- rep(0, length(measure))}
  if (same.pos.allowed == T){same.pos.allowed <- 1} else{
    same.pos.allowed <- 0}
  ac.tot <- rep(NA, length(predict.lat))
  for (i in 1:length(predict.lat)){
    act.lat <- predict.lat[i]; act.lon <- predict.lon[i]
    ind.weight <- 1 / (same.pos.allowed + sqrt((act.lat - lat)^2
                                               + (act.lon - lon)^2))
    ac.tot[i] <- sum(ac_term * ind.weight) / sum(ind.weight)
  }
  return(ac.tot)
}

#function for ac_term
get.1d.ac<-function(resis, ac.sd, lat, long=NULL, contr.fac=NULL, excl.contr.fac=NULL){
  # if no control factor, all is 1
  if(is.null(contr.fac)){contr.fac <- rep(1, length(resis))}
  # if length of excl.contr.fac is null, than every level different
  if(is.null(excl.contr.fac)){excl.contr.fac <- 1:length(resis)}
  # the levels of contr.fac are factors
  contr.fac <- as.factor(contr.fac)
  # the levels of excl.contr.fac are levels
  excl.contr.fac <- as.factor(excl.contr.fac)
  # levels of the factor
  c.fac.lev <- levels(contr.fac)
  # if long is not there, it is set to 0
  if(is.null(long)){long <- rep(0, length(lat))}

  all.ac <- rep(NA, length(resis))
  # for each level of c.fac.lev
  for(i in 1:length(c.fac.lev)){
    # this selects all residuals from the same year,
    # which are not from the respective residual itself
    ind.resis <- resis[contr.fac == c.fac.lev[i] &
                         excl.contr.fac != excl.contr.fac[i]]
    # this selects all latitudes and longitudes from the same year,
    # which are not from the respective residual itself
    ind.lat <- lat[contr.fac == c.fac.lev[i] &
                     excl.contr.fac != excl.contr.fac[i]]
    ind.long <- long[contr.fac == c.fac.lev[i] &
                       excl.contr.fac != excl.contr.fac[i]]
    # what is xx here?
    xx <- unlist(lapply(1:length(ind.resis), function(cell){
      ac.weights = dnorm(x = sqrt((ind.lat[cell] - ind.lat[-cell])^2 +
                                    (ind.long[cell] - ind.long[-cell])^2),
                         mean = 0, sd = ac.sd)
      ac.weights[is.nan(ac.weights)] <- 0
      weighted.mean(ind.resis[-cell], w = ac.weights)
    }))
    all.ac[contr.fac == c.fac.lev[i] & excl.contr.fac != excl.contr.fac[i]] <- xx
  }
  all.ac[is.nan(all.ac)] <- 0
  all.ac[is.na(all.ac)] <- 0
  return(all.ac)
}

#function to
new.get.conf.set <- function(inp.table){#inp.table is expected to comprise two variables: named 'aic' and 'models'
  o <- order(inp.table$aic)
  inp.table <- inp.table[o, ]
  d.aic <- inp.table$aic - min(inp.table$aic)
  e.aic <- exp(-d.aic / 2)
  w.aic <- e.aic / sum(e.aic)
  cum <- cumsum(w.aic)
  c.set <- cum <= 0.95
  c.set[1] <- 1
  m.rank <- rank(inp.table$aic)
  result <- data.frame(inp.table, d.aic, e.aic, w.aic, cum, c.set, m.rank)
  o <- order(o)
  result <- result[o, ]
  #as.data.frame(cbind(o.models,o.aic,d.aic,e.aic, w.aic,cum,c.set,order))
  return(result)
}



#calculate summary stats
#parameter estimates
calculate.mean.coefficients <- function(m_terms, results_res, c_set){

  # !!! ATTENTION
  # here all parameter coeff_, excluding the model
  length_par_est <- length(m_terms)
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

  return(res.par.est.av)
}
