#' @importFrom data.table as.data.table data.table copy melt fread rbindlist setDT setnames setorderv is.data.table ":=" ".N" ".SD" ".BY" "fwrite" "fcase" "between" "uniqueN"
#' @importFrom stats coef fitted predict sd .lm.fit C aggregate cor cov lm mad median na.omit poly qnorm quantile rnorm sigma var
#' @importFrom utils data head tail
#' @importFrom ggplot2 aes aes_string annotate element_blank element_line element_rect element_text expansion facet_grid facet_wrap geom_abline geom_col geom_errorbar geom_histogram geom_hline geom_line geom_point geom_ribbon geom_segment geom_text geom_vline ggplot ggsave ggtitle guide_colorbar guide_legend guides labs scale_alpha_continuous scale_alpha_manual scale_color_discrete scale_color_gradient scale_color_gradient2 scale_color_manual scale_fill_discrete scale_fill_gradient2 scale_fill_manual scale_x_continuous scale_y_continuous theme theme_bw theme_minimal xlab ylab ylim
#' @importFrom grDevices colorRamp rgb
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
#' @title Synthetic Control Methods with Specification Curve Analysis and Machine Learning Interpretability
#' 
#' @description 
#' Researchers using synthetic control methods face a "garden of forking paths" with numerous 
#' specification choices that can dramatically affect results. Existing packages often produce 
#' contradictory estimates even for the same data, creating uncertainty about which analytical 
#' choices drive results. The SCMs package addresses this fundamental challenge through systematic 
#' specification curve analysis integrated with machine learning interpretability methods.
#' 
#' The package provides the first integration of SHAP (SHapley Additive exPlanations) values with 
#' specification curve analysis for synthetic controls. This innovation transforms specification 
#' curves from purely descriptive tools into diagnostic instruments that reveal which analytical 
#' choices drive results. SHAP values quantify how much each specification dimension (constraint 
#' type, outcome model, feature selection) contributes to predicted treatment effects, enabling 
#' researchers to identify the sources of specification sensitivity.
#' 
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Comprehensive Specification Analysis}: Systematic exploration of synthetic control 
#'     modeling choices including constraint types (simplex, lasso, ridge, OLS), outcome models 
#'     (standard, augmented, ridge-augmented), and feature weighting schemes
#'   \item \strong{SHAP Integration}: Machine learning interpretability methods reveal which 
#'     specification choices drive treatment effect estimates, with statistical significance testing
#'   \item \strong{Multiple Estimation Methods}: Supports various synthetic control implementations 
#'     with flexible constraint specifications and optimization approaches  
#'   \item \strong{Advanced Inference}: Implements both Abadie placebo-based inference and bootstrap 
#'     methods for robust statistical testing across specifications
#'   \item \strong{Parallel Processing}: Efficient computation across large specification spaces with 
#'     automatic memory management and progress tracking
#'   \item \strong{Professional Visualization}: Publication-ready specification curve plots with SHAP 
#'     coloring, significance indicators, and customizable aesthetics
#' }
#' 
#' The package workflow proceeds from data preparation through estimation to visualization:
#' data preparation with \code{\link{scdata}}, synthetic control estimation via \code{\link{scest}}, 
#' specification curve generation using \code{\link{spec_curve}}, SHAP analysis with 
#' \code{\link{run_catboost_shap_analysis}}, and visualization through \code{\link{plot_spec_curve}}. 
#' The complete pipeline can be executed with \code{\link{run_spec_curve_analysis}} for turnkey analysis.
#' 
#' \strong{Methodological Innovation:}
#' 
#' Traditional specification curves show the distribution of treatment effect estimates but provide 
#' limited insight into \emph{why} estimates vary. By integrating SHAP values, this package reveals 
#' the causal mechanisms behind specification sensitivity. When a specification produces an unusual 
#' treatment effect, SHAP analysis identifies whether this stems from the constraint type, outcome 
#' model, feature selection, or their interactions.
#' 
#' This diagnostic capability is particularly valuable for synthetic control methods, where 
#' specification choices can interact in complex ways. For example, ridge constraints may perform 
#' well with augmented outcome models but poorly with standard synthetic controls, while simplex 
#' constraints show the opposite pattern. SHAP analysis reveals these interactions automatically.
#' 
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{scdata}}}{Prepare panel data for synthetic control analysis with flexible covariate specifications}
#'   \item{\code{\link{scest}}}{Estimate synthetic control models with various constraint types and optimization methods}
#'   \item{\code{\link{spec_curve}}}{Core computational engine for specification curve analysis across modeling choices}
#'   \item{\code{\link{run_spec_curve_analysis}}}{High-level interface for complete specification curve workflow}
#'   \item{\code{\link{plot_spec_curve}}}{Create publication-ready specification curve visualizations with SHAP integration}
#'   \item{\code{\link{run_catboost_shap_analysis}}}{Compute SHAP values to explain specification choices driving treatment effects}
#' }
#' 
#' @section Example Workflow:
#' \preformatted{
#' # Load and prepare data
#' data <- scdata(df = your_data, 
#'                outcome.var = "gdp",
#'                id.var = "country", 
#'                time.var = "year",
#'                treated.unit = "West Germany",
#'                ...)
#' 
#' # Run complete specification curve analysis  
#' results <- run_spec_curve_analysis(
#'   dataset = your_data,
#'   params = list(
#'     outcome = "gdp",
#'     name_treated_unit = "West Germany",
#'     ...
#'   )
#' )
#' 
#' # Visualize results with SHAP analysis
#' plot_spec_curve(results, 
#'                 name_treated_unit = "West Germany",
#'                 show_shap = TRUE)
#' }
#' 
#' @references
#' Abadie, A., Diamond, A., and Hainmueller, J. (2010). Synthetic Control Methods for Comparative 
#' Case Studies: Estimating the Effect of California's Tobacco Control Program. \emph{Journal of the 
#' American Statistical Association}, 105(490), 493-505.
#' 
#' Ben-Michael, E., Feller, A., and Rothstein, J. (2021). The Augmented Synthetic Control Method. 
#' \emph{Journal of the American Statistical Association}, 116(536), 1789-1803.
#' 
#' Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021). Prediction Intervals for Synthetic Control 
#' Methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.
#' 
#' Lundberg, S. M., and Lee, S. I. (2017). A Unified Approach to Interpreting Model Predictions. 
#' \emph{Advances in Neural Information Processing Systems}, 30, 4765-4774.
#' 
#' @author Zubin Jelveh \email{zjelveh@umd.edu}
#' @keywords synthetic-control specification-curve SHAP machine-learning causal-inference
NULL