#' generate multi-source toy data for regression
#'
#' \code{generateSyntheticData} simualtes random multi-source multi-task regression dataset.
#'
#'
#' @param S is a scaler representing the desired number of sources.
#' @param nY is a vector representing the desired number of samples in each source.
#' @param dX is a scaler representing the desired number of dimensions in X (inputs).
#' @param dY is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.
#' @param ft is a scaler representing the desired number of active features in X.
#' @return A list containing the following elements:
#'  \item{data}{which is a list containing the Y and X data matrices,}
#'  \item{Beta}{the beta parameters used to generate the data.}
#' @export
generateSyntheticData <- function(  S = 2,
                              nY = c(40,80),
                              dX = 50,
                              dY = 1,
                              ft = 10 )
{
  sprintf('    Generating synthetic data from S:%i, dX:%i, dY:%i',S,dX,dY)
  sigma = seq(0.2,0.5,l=S)
  X = matrix(rnorm(sum(nY)*dX),sum(nY),dX)
  betaM = rep(0,dX); betaM[1:ft] = rnorm(ft,0.7,1)

  Beta = matrix(1e-3,S,dX);
  Y = matrix(NA,sum(nY),dY);
  for( s in 1:S){
    Beta[s,1:ft] = abs(rnorm(ft,0.5,0.5)); Beta[s,1:ft] = Beta[s,1:ft]/max(Beta[s,1:ft])
    Beta[s,1:ft] = Beta[s,1:ft]*betaM[1:ft]
    if(s == 1) { st=1; en=nY[1]; }
    else { st=sum(nY[1:(s-1)])+1; en=sum(nY[1:s]); }
    Y[st:en,] =  (X[st:en,] %*% Beta[s,]);
    Y[st:en,] = Y[st:en,] + matrix(sd(Y[st:en,]),nY[s],dY,byrow=T)*sigma[s]*matrix(rnorm(nY[s]*dY,0,1),nY[s],dY)
  }

  data = list(S = S, nY = nY,dX = dX, dY = 1, Y = matrix(Y[,1],nrow(Y),1), X = X, p0 = ft*2)

  return(list(data=data,Beta=Beta))
}

#' demo multi-source model training and interpretative plots
#'
#' \code{demo_bmsr} trains bmsr on random multi-source regression dataset.
#'
#' @export
#' @examples
#'
#' #for multi-source regression
#' demo_bmsr()
#'
demo_bmsr <- function()
{
  file = "bmsr.stan"
  dY = 1
  demo_internal(file = file,dY = dY)
}

#' demo multi-source multi-task model training and interpretative plots
#'
#' \code{demo_bmsmtr} trains bmsmtr on random multi-source multi-task regression dataset.
#'
#' @param dY is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.
#' @export
#' @examples
#'
#' #for multi-source multi-task regression
#' demo_bmsmtr(dY = 3)
#'
demo_bmsmtr <- function(dY = 3)
{
  file = "bmsmtr.stan"
  demo_internal(file = file,dY = dY)
}

#' demo multi-source/multi-task model training and interpretative plots
#'
#' \code{demo_bmsmtr} trains bmsr and bmsmtr on random multi-source/multi-task regression dataset.
#'
#' @param file is the stan file containing the stan code.
#' @param dY is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.
#' @examples
#'
#' #for multi-source multi-task regression
#' demo_internal(file = "bmsr.stan",dY = 1)
#'
#' #for multi-source multi-task regression
#' demo_internal(file = "bmsmtr.stan",dY = 3)
#'
#' @keywords internal
demo_internal <- function(file = "bmsr.stan",dY = 1)
{

  #data
  ll = generateSyntheticData(S = 2,
                             nY = c(40,80),
                             dX = 50,
                             dY = dY,
                             ft = 10);
  data = ll$data; Beta = ll$Beta
  xTest = data$X
  nTest = data$nY
  Y = data$Y
  S = data$S
  dX = ncol(xTest)
  dY = ncol(Y)

  #parameters
  opts = list(iter=1000,seeds=c(12,345,6789),inference="Sampling")
  predFunction = paste0("predict.",file);
  posteriorFunction = paste0("posterior.",file);

  #training
  sprintf('    Training Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
  res = runSTAN(file,data,opts)
  out = res$out

  #prediction
  sprintf('    Predict from Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
  yPred = predictSTAN(predFunction,out,xTest,nTest,yN)
  sprintf('    Performance of Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
  sprintf('    Correlation: %.3f',diag(cor(yPred,Y)))
  sprintf('    RMSE: %.3f',sqrt(mean( (yPred - Y)^2 )) )

  #interpretive plots
  post = getPosterior(posteriorFunction,out)
  tt = c(as.vector(post$beta),as.vector(Beta))
  op <- par(mfrow=c(S,1),mgp=c(0.5,0,0),mar=c(1,2,0.5,0)+0.5)
  cols = c('black','red')
  pchs = c(16,17)
  for(s in 1:S){
    plot(Beta[s,],ylim=c(min(tt),max(tt)),pch=pchs[1],xlab="feature in dX",ylab="weight",tck=0.01,main = paste0('Feature weights in source ',s),col=cols[1],cex.lab=.5, cex.axis=.5, cex.main=.5, cex.sub=.5)
    lines(c(-100,100), c(0,0) ,col="gray"); lines(c(0,0), c(-100,100) ,col="gray")
    points(post$beta[s,],col=cols[2],pch=pchs[2])
    legend("topright",
           legend = c("true", "predicted"),
           col = cols,
           pch = pchs,
           bty = "n",
           cex = 0.7,
           text.col = "black",
           horiz = F)
  }
  par(op)
  sprintf('    COMPLETED Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
}


#' demo real data training and result visualization
#'
#' \code{demo_bmsr_lapatinib_gdsc_ccle} trains bmsr model on real dataset from GDSC and CCLE datasets from the PharmacoGx package.
#' The demo trains a model for lapatinib response predictions.
#'
#' @export
#' @examples
#'
#' #model lapatanib response from GDSC and CCLE
#' demo_bmsr_lapatinib_gdsc_ccle()
#'
demo_bmsr_lapatinib_gdsc_ccle <- function()
{
  model.file = "bmsr.stan";
  #mssr.model = rstan::stan_model(file=model.file)
  if (model.file == 'bmsr.stan')
    mssr.model = stanmodels$bmsr
  if (model.file == 'bmsmtr.stan')
    mssr.model = stanmodels$bmsmtr

  # Download the GDSC and CCLE datasets (using PharmacoGx)
  datasets <- download.and.prepare.datasets()

  # Post-process the expression and drug response data to: (1) exclude all-NA samples, (2) ensure expression and drug response
  # data have the same samples, and (3) standardize both expression and drug response

  expression.matrices <- llply(datasets, .fun = function(dataset) dataset[["expr"]])
  drug.response.matrices <- llply(datasets, .fun = function(dataset) dataset[["response"]])

  # Exclude all-NA samples
  expression.matrices <- llply(expression.matrices, .fun = function(mat) na.omit(mat[, !unlist(lapply(1:ncol(mat), function(i) all(is.na(mat[,i]))))]))
  drug.response.matrices <- llply(drug.response.matrices, .fun = function(mat) mat[, !unlist(lapply(1:ncol(mat), function(i) all(is.na(mat[,i]))))])

  # Ensure samples are consistent between expression and drug response data
  for(ds in names(expression.matrices)) {
    common.samples <- intersect(colnames(expression.matrices[[ds]]), colnames(drug.response.matrices[[ds]]))
    expression.matrices[[ds]] <- expression.matrices[[ds]][, common.samples, drop=F]
    drug.response.matrices[[ds]] <- drug.response.matrices[[ds]][, common.samples, drop=F]
  }

  # Standardize expression and drug response data
  expression.matrices <- llply(expression.matrices, .fun = function(mat) t(scale(t(mat), center = TRUE, scale = TRUE)))
  drug.response.matrices <- llply(drug.response.matrices, .fun = function(mat) t(scale(t(mat), center = TRUE, scale = TRUE)))

  # Extract response for drug to model
  drug <- "Lapatinib"
  drug.response.vectors <- llply(drug.response.matrices, .fun = function(mat) mat[drug, ,drop=FALSE])

  # Perform feature/gene selection using mRMRe (ensemble-based minimum redundancy, maximum relevance)
  # to improve computational efficiency of BMSR
  nms <- names(datasets)
  names(nms) <- nms

  features <-
    llply(nms, .parallel = TRUE,
          .fun = function(ds) {
            new.df <- data.frame(target = as.numeric(drug.response.vectors[[ds]]), t(expression.matrices[[ds]]))
            new.df <- new.df[!is.na(new.df$target), ]
            data <- mRMR.data(data = new.df)
            res <- mRMR.ensemble(data = data, target_indices = 1, feature_count = 50, solution_count = 20)
            genes <- lapply(1:ncol(res@filters[[1]]), function(i) res@feature_names[res@filters[[1]][,i]])
            sort(unique(as.character(unlist(genes))))
          })

  features <- unique(unname(unlist(features)))

  # Limit expression matrices to genes selected by mRMRe
  expression.matrices <- llply(expression.matrices, .fun = function(mat) mat[features, ])

  # Set a random seed
  set.seed(1234)

  # Train BMSR
  fit <- train.mssr(expression.matrices, drug.response.vectors)
  out = fit$out

  # Predict drug response (here using test data = training data)
  ypred <- predict.mssr(expression.matrices, out)

  # Plot predicted versus observed drug response
  prep <- prepare.bmsr.data(expression.matrices, drug.response.vectors)
  yobs <- as.numeric(prep$data.y)
  names(yobs) <- colnames(prep$data.y)

  common.samples <- intersect(names(ypred), names(yobs))
  ypred <- as.numeric(ypred[common.samples])
  yobs <- as.numeric(yobs[common.samples])

  df <- data.frame(ypred = ypred, yobs = yobs)
  g <- ggplot(data = df, aes(x = ypred, y = yobs)) + geom_point() + xlab("Predicted Response") + ylab("Observed Response")

  # Get the posterior means for the beta vectors across both GDSC and CCLE
  posteriorFunction = paste0("posterior.",model.file);
  post = getPosterior(posteriorFunction,out)

  betaShared = as.numeric(post$betaShared)
  names(betaShared) <- features

  # Plot the posterior betas
  df <- data.frame(gene = features, dataset1 = post$beta[1,], dataset2 = post$beta[2,])
  g <- ggplot(data = df, aes(x = dataset1, y = dataset2)) + geom_point()
  g <- g + xlab("GDSC Posterior Beta") + ylab("CCLE Posterior Beta")
}

#' Demo applying BMSR for integrative analysis of Tavor and Beat AML/OHSU datasets to discover biomarkers of venetoclax response.
#'
#' \code{demo_bmsr_venetoclax_tavor_ohsu} trains bmsr model on real dataset from Tavor and Beat AML/OHSU datasets downloaded from
#' ORCESTRA and harmonized with PharmacoGx. The demo trains a model to predict biomarkers of venetoclax response in AML.
#'
#' @export
#' @examples
#'
#' #modelv venetoclax response from Tavor and Beat AML/OHSU datasets
#' demo_bmsr_venetoclax_tavor_ohsu()
#'
demo_bmsr_venetoclax_tavor_ohsu <- function() {

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
  file <- system.file("extdata", file, package = "bmsr")
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
}
