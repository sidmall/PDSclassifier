#' calculateSMI
#'
#' This function calculates the 'Stem Maturation Index' (SMI), a transcriptomic measure that indicates how well a sample or cell aligns along the stem-differentiation scale.
#' It also provides the MYC targets and PRC targets single-sample or single-cell scores.
#'
#' @param data gene expression matrix (with rows as genes) or singlec-cell: seurat/sce objects
#' @param datatype is the data derived from bulk tissue or single-cell (sc)?
#' @param species between 'human' and 'mouse', select one based on whether the input expression data is derived from human or mouse
#'
#' @return a dataframe with rows as samples/cells, MYC Targets, PRC Targets (ssGSEA scores), SMI_unscaled, and SMI.
#' @export
#'
#' @author Sudhir Malla


calculateSMI <- function(data, datatype = c("bulk", "sc"), species = c("human", "mouse")){


  ## Sort out the geneset first ----
  if(species == "mouse") {

    ### Select only mouse geneset ----
    geneset <- myc_prc_gs_list[names(myc_prc_gs_list)[grep("_Mouse", names(myc_prc_gs_list))]]
    names(geneset) <- gsub("_Mouse", "_Targets", names(geneset))

  } else {

    ### Select only human geneset by default ----
    geneset <- myc_prc_gs_list[names(myc_prc_gs_list)[grep("_Human", names(myc_prc_gs_list))]]
    names(geneset) <- gsub("_Human", "_Targets", names(geneset))
  }


  ## get the data: check if it is bulk or single-cell (sc) ----
  if(datatype == "sc") {

    ## NOTE: because the enrichIt function automatically checks for Seurat, SingleCellExperiment or
    ## Matrix, we skip check on datatype

    ## to check, the number of genes matching, we will extract matrix from the single-cell data
    ## using function from escape:::cntEval
    cnts <- cntEval(data)
    ### Check the proportion of genes in each geneset which are present in the expression data
    genes_matched_prop <- lapply(geneset,
                                 function(genes) sum(genes %in% rownames(cnts)) / length(genes))

    ## if less than 5 genes, stop or else calculate scores
    if(all((genes_matched_prop < 0.05) == TRUE)) {
      stop("Less than 5 genes (minimum) from the genesets were found in the data.")
    }
    else {

      ### escape version change
      if(utils::packageVersion("escape") <= "1.6.0") {

        ## run escape to calculate scores ----
        myc_prc_score <- escape::enrichIt(obj = data,
                                          gene.sets = geneset, method = 'ssGSEA',
                                          groups = 10000, cores = 4,
                                          min.size = 5)

      } else {

        ## run escape to calculate scores ----
        myc_prc_score <- escape::escape.matrix(input.data = data,
                                               gene.sets = geneset, method = 'ssGSEA',
                                               groups = 10000,
                                               min.size = 5,
                                               BPPARAM = BiocParallel::MulticoreParam(workers = 4)
                                               )

      }

    }

  }

  if(datatype == "bulk") {

    ## check if the data is in matrix format ----
    if(!is.matrix(data)) {
      stop("data needs to be a matrix with rows as genes.")
    }
    else if(is.matrix(data)) {

      ### Check the proportion of genes in each geneset which are present in the expression data
      genes_matched_prop <- lapply(geneset,
                                   function(genes) sum(genes %in% rownames(data)) / length(genes))

      ## if less than 5 genes, stop or else calculate scores
      if(all((genes_matched_prop < 0.05) == TRUE)) {
        stop("Less than 5 genes (minimum) from the genesets were found in the data.")
      }
      else {

        ### package version GSVA update
        if(utils::packageVersion("GSVA") <= "1.48.3") {

          ## if matrix, calculate ssgsea score ----
          myc_prc_score <- GSVA::gsva(data,
                                      geneset,
                                      min.sz=5,
                                      max.sz=Inf,
                                      verbose = T,
                                      method = 'ssgsea',
                                      ssgsea.norm = F)

        } else{

          ## if matrix, calculate ssgsea score ----
          ssgsea_par <- GSVA::ssgseaParam(data,
                                          geneset,
                                          minSize = 5,
                                          maxSize = Inf,
                                          normalize = F)
          myc_prc_score <- GSVA::gsva(ssgsea_par,
                                      verbose = F)

        }

        ## Take the scores and transpose it ----
        myc_prc_score <- t(myc_prc_score)

      }
    }
  }


  ## into dataframe ----
  myc_prc_score_copy <- as.data.frame(myc_prc_score)

  ## Calculate SMI (unscaled) ----
  myc_prc_score_copy[["SMI_unscaled"]] <- myc_prc_score_copy[,1]-myc_prc_score_copy[,2] ## PRC - MYC
  ## Scale SMI between -1 to 1 ----
  myc_prc_score_copy[["SMI"]] <- scales::rescale(myc_prc_score_copy[["SMI_unscaled"]], to = c(-1, 1))

  ## return myc_prc
  return(myc_prc_score_copy)

}


## The function below was taked from escape R package
cntEval <- function (obj)
{
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- Seurat::GetAssayData(obj)
  }
  else if (inherits(x = obj, what = "SingleCellExperiment")) {
    cnts <- SingleCellExperiment::counts(obj)
  }
  else {
    cnts <- obj
  }
  if (!inherits(cnts, what = "dgCMatrix")) {
    cnts <- Matrix::Matrix(as.matrix(cnts), sparse = TRUE)
  }
  return(cnts)
}
