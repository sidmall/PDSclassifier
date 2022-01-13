#' PDSpredict
#'
#' PDSpredict gives PDS (Pathway Derived Subtypes) calls and prediction scores from input expression data
#'
#' @param x gene expression data in dataframe format with first column as genes
#' @param species between 'human' and 'mouse', select one based on whether the input expression data is derived from human or mouse
#' @param threshold set threshold between 0 to 1 for PDS calls, default 0.6. Those do not meet the threshold are indicated a 'mixed' bioloical subtype.
#'
#' @return a dataframe with Sample ID, PDS prediction scores for each subtye, PDS prediction and PDS prediction with threshold set by user (default 0.6)
#' @export
#'
#' @author Sudhir Malla

PDSpredict <- function(x, species = c("human", "mouse"), threshold = 0.6) {


  ### examine if  '///' string genes are present and any duplicated genes
  if (nrow(x[duplicated(gsub("\\s.*", "", x[, 1])), ]) > 0){
    x <- highMeanGene(x)
  } else {
    rownames(x) <- x[,1]
    x <- x[,-1]
    x <- as.matrix(x)
  }

  ### Get the C2 (BIOCARTA, KEGG, PID, REACTOME) gene set, saved from MSigDB
  if (species == "mouse") {

    c2.geneset <- mouse.c2.geneset.list

  } else {

    c2.geneset <- human.c2.geneset.list
  }

  ### default threshold = 0.6
  if (threshold < 0 | threshold > 1) {
    stop('Set threshold between 0 and 1')
  } else {
    threshold <- threshold}

  ### Generate single sample GSEA scores (ssGSEA)
  message('Calculating ssGSEA scores...')
  ### set seed
  if (as.Date(paste0(R.version$year,'-',R.version$month,'-',R.version$day)) >= as.Date('2019-04-26'))
  {suppressWarnings(set.seed(123, sample.kind = 'Rounding'))} else{set.seed(123)}
  ### ssGSEA
  y <- GSVA::gsva(x, c2.geneset, max.sz = Inf, ## min.size default as it may lead to no gene set scores
                  verbose = F, method = "ssgsea", parallel.sz = 4,
                  ssgsea.norm = F)
  ### Apply scaling with fixed value from training set.
  minmax_diff <- 20928.93
  y.1 <- as.data.frame(sapply(as.data.frame(y), function(s) {s/(minmax_diff)}))
  rownames(y.1) <- rownames(y)
  y.1 <- as.data.frame(y.1)
  y.1[['gs_name']] <- rownames(y.1)
  y.1 <- dplyr::select(y.1, c('gs_name', dplyr::everything()))
  y.1$gs_name <- gsub('\\.', ':', y.1$gs_name)
  cat("ssGSEA Complete!\n")


  ### train data
  train.data.1 <- trainData[,-c(1:2)]
  rownames(train.data.1) <- trainData[,1]
  train.data.1 <- as.data.frame(t(train.data.1))
  train.data.1[['gs_name']] <- rownames(train.data.1)
  train.data.1 <- dplyr::select(train.data.1, c('gs_name', dplyr::everything()))

  ### merge traindata with new data
  merge.1 <- merge(train.data.1, y.1, by = "gs_name")
  rownames(merge.1) <- merge.1[["gs_name"]]
  merge.1 <- merge.1[,-1]
  nrow(train.data.1) == nrow(merge.1)

  ### check for batch effect
  message('Batch Correction (with training set as reference batch)...')
  combi <- suppressMessages(sva::ComBat(dat = as.matrix(merge.1),
                                        batch = c(rep(1,ncol(train.data.1[,-1])),
                                                  rep(2,ncol(y.1[,-1]))),
                                        mod = NULL,
                                        ref.batch = 1,
                                        par.prior = T,
                                        prior.plots = F))

  cat("Batch Correction Complete!\n")

  combi.train.data <- combi[, c(1:ncol(train.data.1[,-1]))]
  new.data.mtx <- combi[,c((ncol(train.data.1[,-1])+1):(ncol(train.data.1[,-1])+ncol(y.1[,-1])))]


  message('Applying the svmRBF classification model...')

  ### Set the seed for reproducibility
  if (as.Date(paste0(R.version$year,'-',R.version$month,'-',R.version$day)) >= as.Date('2019-04-26'))
  {suppressWarnings(set.seed(40218336, sample.kind = 'Rounding'))} else{set.seed(40218336)}

  ### Support Vector Machine with Radial Basis Function Kernal
  ### prediction on new data
  svm.predict <- stats::predict(svm.model, t(new.data.mtx))
  table(svm.predict)

  ### Final decision on prediction calls
  PDS.predict <- stats::predict(svm.model, t(new.data.mtx), type = "prob")
  rownames(PDS.predict) <- rownames(t(new.data.mtx))
  PDS.predict[['Sample_ID']] <- rownames(PDS.predict)
  PDS.predict <- dplyr::select(PDS.predict, c('Sample_ID', dplyr::everything()))
  PDS.predict <- data.frame(PDS.predict, prediction = svm.predict)

  PDS.predict[["PDS_call"]] <- ifelse(PDS.predict[["PDS1"]] > threshold, "PDS1",
                                      ifelse(PDS.predict[["PDS2"]] > threshold, "PDS2",
                                             ifelse(PDS.predict[["PDS3"]] > threshold, "PDS3", "Mixed")))

  cat("Classification Complete!\n")
  beepr::beep(sound = 4)

  PDS.predict <- dplyr::mutate(PDS.predict, PDS_call = factor(PDS.predict$PDS_call,
                                                              levels = c('PDS1', 'PDS2',
                                                                         'PDS3', 'Mixed')))

  print(table(PDS.predict[['PDS_call']]))



  return(PDS.predict)

}
