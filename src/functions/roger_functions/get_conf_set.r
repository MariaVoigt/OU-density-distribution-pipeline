#gets the 95% confidence set from a vector comprising AIC-values
#written by Roger Munrdy
#version from Oct. 29, 2011
#arguments is:
#inp.table: a numerical vector with the aic-values
#returns a data frame comprising the AIC-values handed over (aic),
#the difference between each AIC and the smalles AIC (d.aic)
#Akaike weights (w.aic)
#cumulative Akaike weights (cum)
#whether a model is in the 95% confidence set (c.set; 0 mean no, 1 means yes)
#the rank of the AIC (m.rank)
conf.set<-function(aic){
  o=order(aic)
  #aic=aic[o]
  d.aic=aic-min(aic)
  e.aic=exp(-d.aic/2)
  w.aic= e.aic/sum(e.aic)
  cum=cumsum(w.aic[o])[order(o)]
  c.set=cum<=0.95
  c.set[which.min(aic)]=1
  if (max(cum[c.set==T])<0.95){
    c.set[cum==min(cum[c.set==F])]=T
  }
  m.rank=rank(aic)
  result=data.frame(aic, d.aic, w.aic, cum, c.set, m.rank)
  #o=order(o)
  #result=result[o,]
  return(result)
}
