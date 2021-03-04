suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(SummarizedExperiment))
suppressPackageStartupMessages(p_load(S4Vectors))
suppressPackageStartupMessages(p_load(PharmacoGx))
suppressPackageStartupMessages(p_load(rstan))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(mRMRe))
suppressPackageStartupMessages(p_load(foreach))
suppressPackageStartupMessages(p_load(parallel))
suppressPackageStartupMessages(p_load(ggplot2))

# Integrative modeling of Lapatinib across GDSC and CCLE breast cancer cell lines
# using BMSR for joint analysis over GDSC and CCLE, mRMRe for initial feature selection, and
# PharmacoGx for GDSC and CCLE data access.

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1
## Following is for rstan
options(mc.cores = num.processes)

#source("dataset_utils.R")
#source("bmsr.R")
#source("bmsr_wrappers.R")
require(bmsr)

demo_bmsr_lapatinib_gdsc_ccle()
