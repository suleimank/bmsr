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
suppressPackageStartupMessages(p_load(EnsDb.Hsapiens.v75))
suppressPackageStartupMessages(p_load(edgeR))
suppressPackageStartupMessages(p_load(ggrepel))

# Integrative modeling of Venetoclax across Tavor and OHSU ("Beat AML") AML ex vivo datasets
# using BMSR for joint analysis over Tavor and OHSU, mRMRe for initial feature selection, and
# Orcestra for Tavor and OHSU data acess, and PharmacoGx for data harmonization, drug response
# computation, and expression analysis.

# References:
# "Tavor" dataset: Tavor, Sigal, Tali Shalit, Noa Chapal Ilani, Yoni Moskovitz, Nir Livnat, Yoram Groner, Haim Barr, et al. 2020. “Dasatinib Response in Acute Myeloid Leukemia Is Correlated with FLT3/ITD, PTPN11 Mutations and a Unique Gene Expression Signature.” Haematologica 105 (12): 2795–2804.
# "OHSU" dataset: Tyner, Jeffrey W., Cristina E. Tognon, Daniel Bottomly, Beth Wilmot, Stephen E. Kurtz, Samantha L. Savage, Nicola Long, et al. 2018. “Functional Genomic Landscape of Acute Myeloid Leukaemia.” Nature 562 (7728): 526–31.
# PharmacoGx: Smirnov, Petr, Zhaleh Safikhani, Nehme El-Hachem, Dong Wang, Adrian She, Catharina Olsen, Mark Freeman, et al. 2016. “PharmacoGx: An R Package for Analysis of Large Pharmacogenomic Datasets.” Bioinformatics  32 (8): 1244–46.
# mRMRe: De Jay, Nicolas, Simon Papillon-Cavanagh, Catharina Olsen, Nehme El-Hachem, Gianluca Bontempi, and Benjamin Haibe-Kains. 2013. “mRMRe: An R Package for Parallelized mRMR Ensemble Feature Selection.” Bioinformatics  29 (18): 2365–68.
# ORCESTRA: https://orcestra.ca/
# ORCESTRA: a platform for orchestrating and sharing high-throughput pharmacogenomic analyses (https://www.biorxiv.org/content/10.1101/2020.09.18.303842v1)

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1
## Following is for rstan
options(mc.cores = num.processes)

require(bmsr)

bmsr:::demo_bmsr_venetoclax_tavor_ohsu()

