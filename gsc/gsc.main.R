###############
### gsc main
###############
## preprocess the raw data
source("./gsc/gsc.data.preprocess.R")

# ## if already have the "./gsc/gsc.analysis.cleaned.Rdata", run all the codes below this line
source("./gsc/gsc.data.preparation.R")

## data analysis
source("./gsc/gsc.xmer.jags.R") # Mx, x=1,2,...,20
source("./gsc/gsc.analysis.new.Mass.jags.R") # M21
source("./gsc/gsc.analysis.new.Mass.normal.jags.R") # M22

## posterior summary
source("./gsc/gsc.xmer.post.summary.R") # produce summary statistics