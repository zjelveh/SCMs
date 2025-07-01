#' @importFrom data.table as.data.table data.table copy melt fread rbindlist setDT
#' @importFrom stats coef fitted predict sd
#' @importFrom CVXR Variable Minimize Problem solve
#' @importFrom kernlab alpha
#' @importFrom foreach %dopar%
#' @importFrom ggnewscale new_scale_color
#' @importFrom cowplot plot_grid
#' @importFrom catboost catboost.train catboost.predict catboost.load_pool catboost.get_feature_importance
#' @import methods
#' @import optimx
#' @import stringr
#' @import tidyr
#'
#' @keywords internal
"_PACKAGE"

#' @name SCMs-package
#' @aliases SCMs
#' @title Synthetic Control Methods
#' @description This and that
NULL