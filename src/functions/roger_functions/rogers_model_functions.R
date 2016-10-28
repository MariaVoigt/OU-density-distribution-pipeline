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
get.1d.ac <- function(resis, ac.sd, lat, long = NULL, contr.fac = NULL){
  if (length(contr.fac) == 0){contr.fac = rep(1, length(resis))}
  contr.fac <- as.factor(contr.fac)
  c.fac.lev <- levels(contr.fac)
  if (length(long) == 0){long <- rep(0, length(lat))}
  all.ac <- rep(NA, length(resis))
  for (i in 1:length(c.fac.lev)){
    ind.resis <- resis[contr.fac == c.fac.lev[i]]
    ind.lat <- lat[contr.fac == c.fac.lev[i]]
    ind.long <- long[contr.fac == c.fac.lev[i]]
    xx <- unlist(lapply( 1:length(ind.resis),
                         function(cell){
                           ac.weights <- dnorm(x = sqrt((ind.lat[cell] - ind.lat[-cell])^2 +
                                                          (ind.long[cell] - ind.long[-cell])^2),
                                               mean = 0, sd = ac.sd)
                           #ac.weights is a function
                           ac.weights[is.nan(ac.weights)] = 0
                           weighted.mean(ind.resis[-cell], w = ac.weights)
                         }))
    all.ac[contr.fac == c.fac.lev[i]] <- xx
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
