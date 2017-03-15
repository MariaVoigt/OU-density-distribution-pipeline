#------------------------------------#
# Make bootstrap for abundance model #
#------------------------------------#


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

is_verbose <- options$verbose_script


# input directory
indir <- options$input_directory
if(is_verbose){print(paste("indir", indir))}

# directory in which output is written
outdir <- options$output_directory
if(is_verbose){print(paste("outdir", outdir))}


#---------#
# Globals #
#---------#

indir_fun <- "~/orangutan_density_distribution/src/functions"
if(is_verbose){print(paste("indir_fun", indir_fun))}

cl <- makeForkCluster(outfile = "")
registerDoParallel(cl)

source(file.path(indir_fun, "roger_functions/rogers_model_functions.R"))
source(file.path(indir_fun, "generic/path.to.current.R"))
source(file.path(indir_fun, "roger_functions/get_conf_set.r"))
source(file.path(indir_fun, "roger_functions/helpers.r"))
source(file.path(indir_fun, "roger_functions/diagnostic_fcns.r"))

options("scipen" = 100, "digits" = 4)

#---------------#
#  Import data  #
#---------------#



#load("/homes/mv39zilo/work/Borneo/outreach/Correspondance/November_2016/Roger/images/abundance_model_fitting_2016-12-02.RData")


# include abundMod_results
#oreductirs_obs und m_terms
abundMod_results_path <- path.to.current(indir, "abundMod_results", "rds")
if(is_verbose){print(paste("abundMod_results_path", abundMod_results_path))}
abundMod_results <- readRDS(abundMod_results_path)

predictors_path <- path.to.current(indir, "predictors_observation_scaled", "rds")
if(is_verbose){print(paste("predictors-path", predictors_path))}
predictors_obs <- readRDS(predictors_path)

m_terms_path <- path.to.current(indir, "m_terms", "rds")
if(is_verbose){print(paste("m_terms_path", m_terms_path))}
m_terms <- readRDS(m_terms_path)


ests=apply(abundMod_results [, grepl(x=colnames(abundMod_results ), pattern="coeff")], 2, function(x){
	x[is.na(x)]=0
	sum(x*abundMod_results$w_aic)
})
SEs=apply(abundMod_results [, grepl(x=colnames(abundMod_results ), pattern="SE")], 2, function(x){
	x[is.na(x)]=0
	sum(x*abundMod_results$w_aic)
})
SEs=SEs[-length(SEs)]
names(ests)=gsub(x=names(ests), pattern="coeff_", replacement="", fixed=T)
names(SEs)=gsub(x=names(SEs), pattern="SE_", replacement="", fixed=T)
if(is_verbose){"sum of names of est that is not in SEs"}
sum(names(ests)!=names(SEs))

theta=sum(abundMod_results $theta*abundMod_results$w_aic)
se.theta=sum(abundMod_results $SE.theta*abundMod_results$w_aic)
m.mat=model.matrix(object=as.formula(paste(c("~", paste(c(m_terms, "offset(offset_term)"), collapse="+")), collapse="")), data=predictors_obs)
m.mat=m.mat[, names(ests)]

all.models=paste(abundMod_results $model, "offset(offset_term)", sep="+")
all.models=paste("rv", all.models, sep="~")

if(is_verbose){"finished saving all variables"}

parLapply(cl=cl, X=1:length(cl), function(x){library(MASS)})

n.boots=1000
n.attempts=rep(0, n.boots)
all.boots=matrix(NA, ncol=length(m_terms), nrow=n.boots)
colnames(all.boots)=m_terms
colnames(all.boots)[1]="(Intercept)"

#plot(1, 1, type="n", xlim=c(1, n.boots))
for(i in 1:n.boots){
print(paste("this is the ", i, "boot"))
    xdone=F
	while(!xdone){
		n.attempts[i]=n.attempts[i]+1
		boot.ests=rnorm(n=length(ests), mean=ests, sd=SEs)
		names(boot.ests)=names(ests)
		rv=m.mat%*%boot.ests+predictors_obs$offset_term
		rv=rnbinom(n=nrow(predictors_obs), size=rnorm(n=1, mean=theta, sd=se.theta), mu=exp(rv))
		clusterExport(cl=cl, varlist=c("ests", "SEs", "m.mat", "predictors_obs", "theta", "se.theta", "all.models", "conf.set", "m_terms", "rv"))
		all.mres=parLapply(cl=cl, X=all.models, function(model){
			model <- as.formula(model)
			res <- try(glm.nb(model, data = predictors_obs), silent=T)
			if(class(res)[1]!="try-error"){
				return(list(aic=res$aic, ests=coefficients(res)))
			}else{
				return(list(aic=NA, ests=NA))
			}
		})
		all.aic=unlist(lapply(all.mres, function(x){x$aic}))
		if(sum(is.na(all.aic))==0){
			xdone=T
		}
	}
	w.aic=conf.set(all.aic)$w.aic
	all.ests=matrix(0, ncol=length(m_terms), nrow=length(all.mres))
	colnames(all.ests)=m_terms
	colnames(all.ests)[1]="(Intercept)"
	for(j in 1:length(all.mres)){
		xx=all.mres[[j]]$ests
		all.ests[j, names(xx)]=xx
	}
	xx=apply(all.ests, 2, function(x){sum(x*w.aic)})
	all.boots[i, names(xx)]=xx
#	points(i, 1, col=rainbow(n.boots)[i], pch=19, cex=sqrt(n.attempts[i]))
}


saveRDS(all.boots, file = file.path(outdir, paste0("all_boots_",
                                                   Sys.Date(), ".rds")))

xx=cbind(ests, t(apply(all.boots, 2, quantile, probs=c(0.025, 0.975))))


saveRDS(xx, file = file.path(outdir, paste0("xx_",
                                                   Sys.Date(), ".rds")))


save.image(file.path(outdir, paste0("abundance_model_bootstrap_",
                                    Sys.Date(), ".RData")))

print(paste("Finished model_fitting script, at", Sys.time()))

