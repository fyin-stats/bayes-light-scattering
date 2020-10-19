###################
# lysozyme main
###################
# read and preprocess the raw data
# source("./lysozyme/lysozyme.data.preprocess.R")

# if already have the file "./lysozyme/lysozyme.3d.analysis.cleaned.Rdata", run all the codes below this line
source("./lysozyme/lysozyme.data.preparation.R")

# data analysis
source("./lysozyme/lysozyme.complete.pool.constant.dndc.noconmeas.R")
source("./lysozyme/lysozyme.complete.pool.constant.dndc.noconmeas.diffuse.R")
source("./lysozyme/lysozyme.complete.pool.constant.dndc.noconmeas.dimer.R")
source("./lysozyme/lysozyme.complete.pool.constant.dndc.noconmeas.M12.R")
source("./lysozyme/lysozyme.complete.pool.constant.dndc.noconmeas.Mass.R")
source("./lysozyme/lysozyme.complete.pool.constant.dndc.noconmeas.no.adjustment.R")
# figures, summary statistics
source("./lysozyme/lysozyme.berkman.post.summary.R")