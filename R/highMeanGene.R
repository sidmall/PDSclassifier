#' highMeanGene
#'
#' Checks for duplicated genes from input expression data, selects gene with highest mean value from
#' duplicated genes, and produces expression matrix with uniques genes as rows.
#'
#' @param df dataframe with first column as genes
#'
#' @return expression matrix with unique genes as rows
#' @export

highMeanGene <- function (df) {
  x.df <- as.data.frame(df)
  x.df <- dplyr::mutate(x.df, gene_noString = gsub("\\s.*",
                                                   "", x.df[, 1]))
  x.df <- dplyr::select(x.df, "gene_noString", dplyr::everything())
  x.df <- dplyr::mutate(x.df, mean = apply(x.df[, -c(1:2)],
                                           1, mean))
  x.df <- dplyr::select(x.df, "mean", dplyr::everything())
  x.df <- x.df[order(x.df[["mean"]], decreasing = T), ]
  x.df <- x.df[-which(duplicated(x.df[["gene_noString"]])),
               ]
  x.mtx <- as.matrix(x.df[, -c(1:3)])
  rownames(x.mtx) <- x.df[["gene_noString"]]
  x.mtx <- x.mtx[order(rownames(x.mtx), decreasing = F), ]
  return(x.mtx)
}
