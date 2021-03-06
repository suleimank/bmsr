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

model.file = "bmsr.stan";

# Download the PSets for OHSU/Beat AML and Tavor from ORCESTRA
# We will then access these PSets through PharmacoGx functions
orcestra.beat.aml.pset.link <- "https://zenodo.org/record/4582786/files/BeatAML.rds?download=1"
orcestra.tavor.pset.link <- "https://zenodo.org/record/4585705/files/Tavor.rds?download=1"

pset.links <- list("OHSU" = orcestra.beat.aml.pset.link, "Tavor" = orcestra.tavor.pset.link)

psets <- download.psets(pset.links)

datasets <- names(psets)
names(datasets) <- datasets

## Extract expression data from PharmacoGx
expression.matrices <-
  llply(psets,
        .fun = function(pset) {
                 molecularProfiles(pset, mDataType="rnaseq")
	       })

## OSHU data have Ensembl ids; let's translate to symbols
translation.tbl <- data.frame(ensg.full = rownames(expression.matrices[["OHSU"]]),
                              ensg = unlist(lapply(rownames(expression.matrices[["OHSU"]]), function(str) gsub(str, pattern="\\.(\\d)+$", replacement=""))))

gene_info <- ensembldb::genes(EnsDb.Hsapiens.v75, filter=GeneIdFilter(translation.tbl$ensg), return.type='data.frame') %>%
                dplyr::select(gene_id, gene_name, gene_biotype)  

translation.tbl <- merge(translation.tbl, gene_info, by.x = c("ensg"), by.y = c("gene_id"), all = FALSE)

# Let's limit to genes that are highly expressed and highly variable in AML
file <- "aml-highly-expressed-variable-genes.tsv"
aml.genes <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)$gene
translation.tbl <- subset(translation.tbl, ensg %in% aml.genes)

# Remove any duplicate mappings (lazy!)
my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
translation.tbl <- subset(translation.tbl, !my.dup(ensg.full))
translation.tbl <- subset(translation.tbl, !my.dup(ensg))
translation.tbl <- subset(translation.tbl, !my.dup(gene_name))

# Translate OHSU gene ids to symbols
common.ensg.genes <- intersect(rownames(expression.matrices[["OHSU"]]), translation.tbl$ensg.full)
translation.tbl <- subset(translation.tbl, ensg.full %in% common.ensg.genes)
expression.matrices[["OHSU"]] <- expression.matrices[["OHSU"]][rownames(expression.matrices[["OHSU"]]) %in% common.ensg.genes,]
rownames(translation.tbl) <- translation.tbl$ensg.full
rownames(expression.matrices[["OHSU"]]) <- translation.tbl[rownames(expression.matrices[["OHSU"]]), "gene_name"]

# Translate Tavor counts to log cpms
expression.matrices[["Tavor"]] <- cpm(expression.matrices[["Tavor"]], log = TRUE)

# Extract area above the curve drug response from PharmacoGx
drug.response.matrices <-
  llply(psets,
        .fun = function(pset) {
	         summarizeSensitivityProfiles(pset,
                                              sensitivity.measure='aac_recomputed',
                                              summary.stat="median",
                                              verbose=FALSE)
               })


# Exclude all-NA samples
expression.matrices <- llply(expression.matrices, .fun = function(mat) na.omit(mat[, !unlist(lapply(1:ncol(mat), function(i) all(is.na(mat[,i]))))]))
drug.response.matrices <- llply(drug.response.matrices, .fun = function(mat) mat[, !unlist(lapply(1:ncol(mat), function(i) all(is.na(mat[,i]))))])

# Ensure samples are consistent between expression and drug response data
for(ds in names(expression.matrices)) {
  common.samples <- intersect(colnames(expression.matrices[[ds]]), colnames(drug.response.matrices[[ds]]))
  expression.matrices[[ds]] <- expression.matrices[[ds]][, common.samples, drop=F]
  drug.response.matrices[[ds]] <- drug.response.matrices[[ds]][, common.samples, drop=F]
}

# Limit to common genes across datasets
common.genes <- Reduce(intersect, lapply(expression.matrices, function(mat) rownames(mat)))
expression.matrices <- llply(expression.matrices, .fun = function(mat) mat[common.genes,])

# Standardize expression and drug response data
expression.matrices <- llply(expression.matrices, .fun = function(mat) t(scale(t(mat), center = TRUE, scale = TRUE)))
drug.response.matrices <- llply(drug.response.matrices, .fun = function(mat) t(scale(t(mat), center = TRUE, scale = TRUE)))

# Extract response for drug to model
drug <- "Venetoclax"
drug.response.vectors <- llply(drug.response.matrices, .fun = function(mat) mat[drug, ,drop=FALSE])

perform.feature.selection <- TRUE

if(perform.feature.selection) {

  # Set a random seed
  set.seed(1234)
  
  # Perform feature/gene selection using mRMRe (ensemble-based minimum redundancy, maximum relevance)
  # to improve computational efficiency of BMSR
  features <-
    llply(datasets, .parallel = TRUE,
          .fun = function(ds) {
                   new.df <- data.frame(target = as.numeric(drug.response.vectors[[ds]]), t(expression.matrices[[ds]]))
                   new.df <- new.df[!is.na(new.df$target), ]
                   data <- mRMR.data(data = new.df)
                   # res <- mRMR.ensemble(data = data, target_indices = 1, feature_count = 200, solution_count = 20)
  		   res <- mRMR.classic(data = data, target_indices = 1, feature_count = 100)
                   genes <- lapply(1:ncol(res@filters[[1]]), function(i) res@feature_names[res@filters[[1]][,i]])
                   genes <- sort(unique(as.character(unlist(genes))))
  		   genes
                 })

  all.features <- unique(unname(unlist(features)))
  common.genes <- Reduce(intersect, lapply(expression.matrices, function(mat) rownames(mat)))
  all.features <- intersect(all.features, common.genes)

  # Limit expression matrices to genes selected by mRMRe
  # expression.matrices <- llply(expression.matrices, .fun = function(mat) mat[rownames(mat) %in% all.features, ])
  expression.matrices <- llply(expression.matrices, .fun = function(mat) mat[all.features,])
}

# Set a random seed
set.seed(1234)

# Train BMSR
fit <- train.mssr(expression.matrices, drug.response.vectors)
out = fit$out

# Get the posterior means for the beta vectors across both Tavor and OHSU
posteriorFunction = paste0("posterior.",model.file);
post = getPosterior(posteriorFunction,out)

betaShared = as.numeric(post$betaShared)
names(betaShared) <- features

# Plot the posterior betas
df <- data.frame(gene = all.features, t(post$beta))
colnames(df)[2:(length(datasets)+1)] <- names(datasets)

# Manually select a few outliers
df.outlier <- subset(df, Tavor < -0.035 | Tavor > 0.05)

# o <- order(df$dataset2, decreasing=FALSE)
# df <- df[o, ]

g <- ggplot(data = df, aes(x = OHSU, y = Tavor)) + geom_point()
g <- g + xlab("OHSU Posterior Beta") + ylab("Tavor Posterior Beta")
g <- g + theme(text = element_text(size = 15))
g <- g + geom_text_repel(data = df.outlier, aes(x = OHSU, y = Tavor, label = gene))

g
ggsave("demo_bmsr_venetoclax_tavor_ohsu.pdf")
ggsave("demo_bmsr_venetoclax_tavor_ohsu.png")
