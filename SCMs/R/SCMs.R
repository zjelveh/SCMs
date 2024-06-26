###############################################################################

#' @title Prediction Intervals for Synthetic Control Methods
#'
#' @description The command implements estimation and inference procedures for Synthetic Control (SC) methods using least squares, lasso, ridge, or simplex-type constraints. Uncertainty is quantified using prediction
#' intervals according to \href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, Feng, and Titiunik (2021)}. \code{\link{scpi}} returns the estimated
#' post-treatment series for the synthetic unit through the command \code{\link{scest}} and quantifies in-sample and out-of-sample uncertainty to provide confidence intervals
#' for each point estimate.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are:  \link{scdata} and \link{scdataMulti} for data preparation in the single and multiple treated unit(s) cases, respectively,
#' \link{scest} for point estimation, \link{scplot} and \link{scplotMulti} for plots in the single and multiple treated unit(s) cases, respectively.
#'
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#'
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#'
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#'
#' @param data a class 'scdata' object, obtained by calling \code{\link{scdata}}, or class 'scdataMulti' obtained via \code{\link{scdataMulti}}.
#' @param w.constr a list specifying the constraint set the estimated weights of the donors must belong to.
#' \code{w.constr} can contain up to five elements:
#' - `\code{p}', a scalar indicating the norm to be used (\code{p} should be one of "no norm", "L1", and "L2")
#' - `\code{dir}', a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
#' - `\code{Q}', a scalar defining the value of the constraint on the norm
#' - `\code{lb}', a scalar defining the lower bound on the weights. It can be either 0 or \code{-Inf}.
#' - `\code{name}', a character selecting one of the default proposals
#' See the \strong{Details} section for more.
#' @param P a \eqn{I\cdot T_1\times I\cdot (J+KM)} matrix containing the design matrix to be used to obtain the predicted.
#' post-intervention outcome of the synthetic control unit. \eqn{T_1} is the number of post-treatment periods,
#' \eqn{J} is the size of the donor pool, and \eqn{K_1} is the number of covariates used for adjustment in the outcome equation.
#' @param V specifies the type of weighting matrix to be used when minimizing the sum of squared residuals
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w})}
#' The default is the identity matrix, so equal weight is given to all observations. In the case of multiple treated observations
#' (you used \code{\link{scdataMulti}} to prepare the data), the user can specify \code{V} as a string equal to either "separate" or "pooled".
#' If \code{scdata()} was used to prepare the data, \code{V} is automatically set to "separate" as the two options are
#' equivalent. See the \strong{Details} section for more.
#' @param V.mat A conformable weighting matrix \eqn{\mathbf{V}} to be used in the minimization of the sum of squared residuals
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r}).}
#' See the \strong{Details} section for more information on how to prepare this matrix.
#' @param solver a string containing the name of the solver used by \code{CVXR} when computing the weights. You can check which solvers are available
#' on your machine by running \code{CVXR::installed_solvers()}. More information on what different solvers do can be found
#' at the following link https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/. "OSQP" is the default solver when 'lasso'
#' is the constraint type, whilst "ECOS" is the default in all other cases.
#' @param rho a string specifying the regularizing parameter that imposes sparsity on the estimated vector of weights. If
#' \code{rho = 'type-1'} (the default), then the tuning parameter is computed based on optimization inequalities. Users can provide a scalar
#' with their own value for \code{rho}. Other options are described in the \strong{Details} section.
#' @param rho.max a scalar indicating the maximum value attainable by the tuning parameter \code{rho}.
#' @param lgapp selects the way local geometry is approximated in simulation. The options are "generalized"
#' and "linear". The first one accommodates for possibly non-linear constraints, whilst the second one is valid
#' with linear constraints only.
#' @param u.missp a logical indicating if misspecification should be taken into account when dealing with \eqn{\mathbf{u}}.
#' @param u.order a scalar that sets the order of the polynomial in \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{u}}.
#' The default is \code{u.order = 1}, however if there is risk of over-fitting, the command automatically sets it
#' to \code{u.order = 0}. See the \strong{Details} section for more information.
#' @param u.lags a scalar that sets the number of lags of \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{u}}.
#' The default is \code{u.lags = 0}, however if there is risk of over-fitting, the command automatically sets it
#' to \code{u.lags = 0}. See the \strong{Details} section for more information.
#' @param u.design a matrix with the same number of rows of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} and whose columns specify the design matrix
#' to be used when modeling the estimated pseudo-true residuals \eqn{\mathbf{u}}.
#' @param u.sigma a string specifying the type of variance-covariance estimator to be used when estimating
#'    the conditional variance of \eqn{\mathbf{u}}.
#' @param u.alpha a scalar specifying the confidence level for in-sample uncertainty, i.e. 1 - \code{u.alpha} is the confidence level.
#' @param e.method a string selecting the method to be used in quantifying out-of-sample uncertainty among:
#'  "gaussian" which uses conditional subgaussian bounds; "ls" which specifies a location-scale model for \eqn{\mathbf{u}}; "qreg" which employs a
#'  quantile regressions to get the conditional bounds; "all" uses each one of the previous methods.
#' @param e.order a scalar that sets the order of the polynomial in \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{e}}.
#' The default is \code{e.order = 1}, however if there is risk of over-fitting, the command automatically sets it
#' to \code{e.order = 0}. See the \strong{Details} section for more information.
#' @param e.lags a scalar that sets the number of lags of \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{e}}.
#' The default is \code{e.order = 1}, however if there is risk of over-fitting, the command automatically sets it
#' to \code{e.order = 0}. See the \strong{Details} section for more information.
#' @param e.design a matrix with the same number of rows of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} and whose columns specify the design matrix
#' to be used when modeling the estimated out-of-sample residuals \eqn{\mathbf{e}}.
#' @param e.alpha a scalar specifying the confidence level for out-of-sample uncertainty, i.e. 1 - \code{e.alpha} is the confidence level.
#' @param sims a scalar providing the number of simulations to be used in quantifying in-sample uncertainty.
#' @param plot a logical specifying whether \code{\link{scplot}} should be called and a plot saved in the current working
#' directory. For more options see \code{\link{scplot}}.
#' @param plot.name a string containing the name of the plot (the format is by default .png). For more options see \code{\link{scplot}}.
#' @param cores number of cores to be used by the command. The default is one.
#' @param w.bounds a \eqn{N_1\cdot T_1\times 2} matrix with the user-provided bounds on \eqn{\beta}. If \code{w.bounds} is provided, then
#' the quantification of in-sample uncertainty is skipped. It is possible to provide only the lower bound or the upper bound
#' by filling the other column with \code{NA}s.
#' @param e.bounds a \eqn{N_1\cdot T_1\times 2} matrix with the user-provided bounds on \eqn{(\widehat{\mathbf{w}},
#' \widehat{\mathbf{r}})^{\prime}}. If \code{e.bounds} is provided, then
#' the quantification of out-of-sample uncertainty is skipped. It is possible to provide only the lower bound or the upper bound
#' by filling the other column with \code{NA}s.
#'
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used to produce the plot.
#'
#' @param verbose if \code{TRUE} prints additional information in the console.
#'
#' @return
#' The function returns an object of class 'scpi' containing three lists. The first list is labeled 'data' and contains used
#' data as returned by \code{\link{scdata}} and some other values.
#' \item{A}{a matrix containing pre-treatment features of the treated unit(s).}
#' \item{B}{a matrix containing pre-treatment features of the control units.}
#' \item{C}{a matrix containing covariates for adjustment.}
#' \item{P}{a matrix whose rows are the vectors used to predict the out-of-sample series for the synthetic unit(s).}
#' \item{Y.pre}{a matrix containing the pre-treatment outcome of the treated unit(s).}
#' \item{Y.post}{a matrix containing the post-treatment outcome of the treated unit(s).}
#' \item{Y.pre.agg}{a matrix containing the aggregate pre-treatment outcome of the treated unit(s). This differs from
#' Y.pre only in the case 'effect' in \code{scdataMulti()} is set to either 'unit' or 'time'.}
#' \item{Y.post.agg}{a matrix containing the aggregate post-treatment outcome of the treated unit(s). This differs from
#' Y.post only in the case 'effect' in \code{scdataMulti()} is set to either 'unit' or 'time'.}
#' \item{Y.donors}{a matrix containing the pre-treatment outcome of the control units.}
#' \item{specs}{a list containing some specifics of the data:
#' \itemize{
#' \item{\code{J}, the number of control units}
#' \item{\code{K}, a numeric vector with the number of covariates used for adjustment for each feature}
#' \item{\code{M}, number of features}
#' \item{\code{KM}, the total number of covariates used for adjustment}
#' \item{\code{KMI}, the total number of covariates used for adjustment}
#' \item{\code{I}, number of treated unit(s)}
#' \item{\code{period.pre}, a numeric vector with the pre-treatment period}
#' \item{\code{period.post}, a numeric vector with the post-treatment period}
#' \item{\code{T0.features}, a numeric vector with the number of periods used in estimation for each feature}
#' \item{\code{T1.outcome}, the number of post-treatment periods}
#' \item{\code{constant}, for internal use only}
#' \item{\code{effect}, for internal use only}
#' \item{\code{anticipation}, number of periods of potential anticipation effects}
#' \item{\code{out.in.features}, for internal use only}
#' \item{\code{treated.units}, list containing the IDs of all treated units}
#' \item{\code{donors.list}, list containing the IDs of the donors of each treated unit}}}
#'
#' The second list is labeled 'est.results' containing all the results from \code{\link{scest}}.
#' \item{w}{a matrix containing the estimated weights of the donors.}
#' \item{r}{a matrix containing the values of the covariates used for adjustment.}
#' \item{b}{a matrix containing \eqn{\mathbf{w}} and \eqn{\mathbf{r}}.}
#' \item{Y.pre.fit}{a matrix containing the estimated pre-treatment outcome of the SC unit(s).}
#' \item{Y.post.fit}{a matrix containing the estimated post-treatment outcome of the SC unit(s).}
#' \item{A.hat}{a matrix containing the predicted values of the features of the treated unit(s).}
#' \item{res}{a matrix containing the residuals \eqn{\mathbf{A}-\widehat{\mathbf{A}}}.}
#' \item{V}{a matrix containing the weighting matrix used in estimation.}
#' \item{w.constr}{a list containing the specifics of the constraint set used on the weights.}
#'
#' The third list is labeled 'inference.results' and contains all the inference-related results.
#' \item{CI.in.sample}{a matrix containing the prediction intervals taking only in-sample uncertainty in to account.}
#' \item{CI.all.gaussian}{a matrix containing the prediction intervals estimating out-of-sample uncertainty with sub-Gaussian bounds.}
#' \item{CI.all.ls}{a matrix containing the prediction intervals estimating out-of-sample uncertainty with a location-scale model.}
#' \item{CI.all.qreg}{a matrix containing the prediction intervals estimating out-of-sample uncertainty with quantile regressions.}
#' \item{bounds}{a list containing the estimated bounds (in-sample and out-of-sample uncertainty).}
#' \item{Sigma}{a matrix containing the estimated (conditional) variance-covariance \eqn{\boldsymbol{\Sigma}}.}
#' \item{u.mean}{a matrix containing the estimated (conditional) mean of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{u.var}{a matrix containing the estimated (conditional) variance-covariance of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{e.mean}{a matrix containing the estimated (conditional) mean of the out-of-sample error \eqn{e}.}
#' \item{e.var}{a matrix containing the estimated (conditional) variance of the out-of-sample error \eqn{e}.}
#' \item{u.missp}{a logical indicating whether the model has been treated as misspecified or not.}
#' \item{u.lags}{an integer containing the number of lags in B used in predicting moments of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{u.order}{an integer containing the order of the polynomial in B used in predicting moments of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{u.sigma}{a string indicating the estimator used for \code{Sigma}.}
#' \item{u.user}{a logical indicating whether the design matrix to predict moments of \eqn{\mathbf{u}} was user-provided.}
#' \item{u.T}{a scalar indicating the number of observations used to predict moments of \eqn{\mathbf{u}}.}
#' \item{u.params}{a scalar indicating the number of parameters used to predict moments of \eqn{\mathbf{u}}.}
#' \item{u.D}{the design matrix used to predict moments of \eqn{\mathbf{u}},}
#' \item{u.alpha}{a scalar determining the confidence level used for in-sample uncertainty, i.e. 1-\code{u.alpha} is the confidence level.}
#' \item{e.method}{a string indicating the specification used to predict moments of the out-of-sample error \eqn{e}.}
#' \item{e.lags}{an integer containing the number of lags in B used in predicting moments of the out-of-sample error \eqn{e}.}
#' \item{e.order}{an integer containing the order of the polynomial in B used in predicting moments of the out-of-sample error \eqn{e}.}
#' \item{e.user}{a logical indicating whether the design matrix to predict moments of \eqn{e} was user-provided.}
#' \item{e.T}{a scalar indicating the number of observations used to predict moments of \eqn{\mathbf{u}}.}
#' \item{e.params}{a scalar indicating the number of parameters used to predict moments of \eqn{\mathbf{u}}.}
#' \item{e.alpha}{a scalar determining the confidence level used for out-of-sample uncertainty, i.e. 1-\code{e.alpha} is the confidence level.}
#' \item{e.D}{the design matrix used to predict moments of \eqn{\mathbf{u}},}
#' \item{rho}{an integer specifying the estimated regularizing parameter that imposes sparsity on the estimated vector of weights.}
#' \item{Q.star}{a list containing the regularized constraint on the norm.}
#' \item{epskappa}{a vector containing the estimates for \eqn{\epsilon_{\kappa}}.}
#' \item{sims}{an integer indicating the number of simulations used in quantifying in-sample uncertainty.}
#' \item{failed.sims}{a matrix containing the number of failed simulations per post-treatment period to estimate lower and upper bounds.}
#'
#'
#' @details
#' Information is provided for the simple case in which \eqn{N_1=1} if not specified otherwise.
#' \itemize{
#' \item{\strong{Estimation of Weights.} \code{w.constr} specifies the constraint set on the weights. First, the element
#' \code{p} allows the user to choose between imposing a constraint on either the L1 (\code{p = "L1"})
#' or the L2 (\code{p = "L2"}) norm of the weights and imposing no constraint on the norm (\code{p = "no norm"}).
#' Second, \code{Q} specifies the value of the constraint on the norm of the weights.
#' Third, \code{lb} sets the lower bound of each component of the vector of weights.
#' Fourth, \code{dir} sets the direction of the constraint on the norm in case \code{p = "L1"}
#' or \code{p = "L2"}. If \code{dir = "=="}, then
#' \deqn{||\mathbf{w}||_p = Q,\:\:\: w_j \geq lb,\:\: j =1,\ldots,J}
#' If instead \code{dir = "<="}, then
#' \deqn{||\mathbf{w}||_p \leq Q,\:\:\: w_j \geq lb,\:\: j =1,\ldots,J}
#' If instead \code{dir = "NULL"} no constraint on the norm of the weights is imposed.
#'
#' An alternative to specifying an ad-hoc constraint set on the weights would be
#' choosing among some popular types of constraints. This can be done by including the element
#' `\code{name}' in the list \code{w.constr}. The following are available options:
#' \itemize{
#' \item {If \code{name == "simplex"} (the default), then
#' \deqn{||\mathbf{w}||_1 = 1,\:\:\: w_j \geq 0,\:\: j =1,\ldots,J.}}
#'
#' \item {If \code{name == "lasso"}, then
#' \deqn{||\mathbf{w}||_1 \leq Q,}
#' where \code{Q} is by default equal to 1 but it can be provided as an element of the list (eg. \code{w.constr =
#' list(name = "lasso", Q = 2)}).}
#'
#' \item{If \code{name == "ridge"}, then
#' \deqn{||\mathbf{w}||_2 \leq Q,}
#' where \eqn{Q} is a tuning parameter that is by default computed as
#' \deqn{(J+KM) \widehat{\sigma}_u^{2}/||\widehat{\mathbf{w}}_{OLS}||_{2}^{2}}
#' where \eqn{J} is the number of donors and \eqn{KM} is the total number of covariates used for adjustment.
#' The user can provide \code{Q} as an element of the list (eg. \code{w.constr =
#' list(name = "ridge", Q = 1)}).}
#'
#' \item{If \code{name == "ols"}, then the problem is unconstrained and the vector of weights
#' is estimated via ordinary least squares.}
#'
#' \item{If \code{name == "L1-L2"}, then
#' \deqn{||\mathbf{w}||_1 = 1,\:\:\: ||\mathbf{w}||_2 \leq Q,}
#' where \eqn{Q} is a tuning parameter computed as in the "ridge" case.}
#' }}
#'
#' \item{\strong{Weighting Matrix.}
#' \itemize{
#' \item{if \code{V <- "separate"}, then \eqn{\mathbf{V} = \mathbf{I}} and the minimized objective function is
#' \deqn{\sum_{i=1}^{N_1} \sum_{l=1}^{M} \sum_{t=1}^{T_{0}}\left(a_{t, l}^{i}-\mathbf{b}_{t, l}^{{i \prime }} \mathbf{w}^{i}-\mathbf{c}_{t, l}^{{i \prime}} \mathbf{r}_{l}^{i}\right)^{2},}
#' which optimizes the separate fit for each treated unit.}
#' \item{if \code{V <- "pooled"}, then \eqn{\mathbf{V} = \mathbf{1}\mathbf{1}'\otimes \mathbf{I}} and the minimized objective function is
#' \deqn{\sum_{l=1}^{M} \sum_{t=1}^{T_{0}}\left(\frac{1}{N_1^2} \sum_{i=1}^{N_1}\left(a_{t, l}^{i}-\mathbf{b}_{t, l}^{i \prime} \mathbf{w}^{i}-\mathbf{c}_{t, l}^{i\prime} \mathbf{r}_{l}^{i}\right)\right)^{2},}
#' which optimizes the pooled fit for the average of the treated units.}
#' \item{if the user wants to provide their own weighting matrix, then it must use the option \code{V.mat} to input a \eqn{v\times v} positive-definite matrix, where \eqn{v} is the
#' number of rows of \eqn{\mathbf{B}} (or \eqn{\mathbf{C}}) after potential missing values have been removed. In case the user
#' wants to provide their own \code{V}, we suggest to check the appropriate dimension \eqn{v} by inspecting the output
#' of either \code{scdata} or \code{scdataMulti} and check the dimensions of \eqn{\mathbf{B}} (and \eqn{\mathbf{C}}). Note that
#' the weighting matrix could cause problems to the optimizer if not properly scaled. For example, if \eqn{\mathbf{V}} is diagonal
#' we suggest to divide each of its entries by \eqn{\|\mathrm{diag}(\mathbf{V})\|_1}.}
#' }}
#'
#'
#' \item{\strong{Regularization.} \code{rho} is estimated through the formula
#' \deqn{\varrho = \mathcal{C}\frac{\log (T_0)^c}{T_0^{1/2}}}
#' where \eqn{\mathcal{C} = \widehat{\sigma}_u / \min_j \widehat{\sigma}_{b_j}} if \code{rho = 'type-1'},
#' \eqn{\mathcal{C} = \max_{j}\widehat{\sigma}_{b_j}\widehat{\sigma}_{u} / \min_j \widehat{\sigma}_{b_j}^2} if \code{rho = 'type-2'}, and
#' \eqn{\mathcal{C} = \max_{j}\widehat{\sigma}_{b_ju} / \min_j \widehat{\sigma}_{b_j}^2} if \code{rho = 'type-3'},
#'
#' \code{rho} defines a new sparse weight vector as
#'  \deqn{\widehat{w}^\star_j = \mathbf{1}(\widehat{w}_j\geq \varrho)}
#' }
#'
#' \item{\strong{In-sample uncertainty.} To quantify in-sample uncertainty it is necessary to model the pseudo-residuals \eqn{\mathbf{u}}.
#' First of all, estimation of the first moment of \eqn{\mathbf{u}} can be controlled through
#' the option \code{u.missp}. When \code{u.missp = FALSE}, then \eqn{\mathbf{E}[u\: |\: \mathbf{D}_u]=0}. If instead \code{u.missp = TRUE},
#' then \eqn{\mathbf{E}[\mathbf{u}\: |\: \mathbf{D}_u]} is estimated using a linear regression of
#' \eqn{\widehat{\mathbf{u}}} on \eqn{\mathbf{D}_u}. The default set of variables in \eqn{\mathbf{D}_u} is composed of \eqn{\mathbf{B}},
#' \eqn{\mathbf{C}} and, if required, it is augmented with lags (\code{u.lags}) and polynomials (\code{u.order}) of \eqn{\mathbf{B}}.
#' The option \code{u.design} allows the user to provide an ad-hoc set of variables to form \eqn{\mathbf{D}_u}.
#' Regarding the second moment of \eqn{\mathbf{u}}, different estimators can be chosen:
#' HC0, HC1, HC2, HC3, and HC4 using the option \code{u.sigma}.}
#'
#' \item{\strong{Out-of-sample uncertainty.} To quantify out-of-sample uncertainty it is necessary to model the out-of-sample residuals
#' \eqn{\mathbf{e}} and estimate relevant moments. By default, the design matrix used during estimation \eqn{\mathbf{D}_e} is composed of the blocks in
#' \eqn{\mathbf{B}} and \eqn{\mathbf{C}} corresponding to the outcome variable. Moreover, if required by the user, \eqn{\mathbf{D}_e}
#' is augmented with lags (\code{e.lags}) and polynomials (\code{e.order}) of \eqn{\mathbf{B}}. The option \code{e.design} allows the user to provide an
#' ad-hoc set of variables to form \eqn{\mathbf{D}_e}. Finally, the option \code{e.method} allows the user to select one of three
#' estimation methods: "gaussian" relies on conditional sub-Gaussian bounds; "ls" estimates conditional bounds using a location-scale
#' model; "qreg" uses conditional quantile regression of the residuals \eqn{\mathbf{e}} on \eqn{\mathbf{D}_e}.}
#'
#' \item{\strong{Residual Estimation Over-fitting.} To estimate conditional moments of \eqn{\mathbf{u}} and \eqn{e_t}
#' we rely on two design matrices, \eqn{\mathbf{D}_u} and \eqn{\mathbf{D}_e} (see above). Let \eqn{d_u} and \eqn{d_e} be the number of
#' columns in \eqn{\mathbf{D}_u} and \eqn{\mathbf{D}_e}, respectively. Assuming no missing values and balanced features, the
#' number of observation used to estimate moments of \eqn{\mathbf{u}} is \eqn{N_1\cdot T_0\cdot M}, whilst for moments of \eqn{e_t} is \eqn{T_0}.
#' Our rule of thumb to avoid over-fitting is to check if \eqn{N_1\cdot T_0\cdot M \geq d_u + 10} or \eqn{T_0 \geq d_e + 10}. If the
#' former condition is not satisfied we automatically set \code{u.order = u.lags = 0}, if instead the latter is not met
#' we automatically set \code{e.order = e.lags = 0}.}
#'
#' }
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Yingjie Feng, Tsinghua University. \email{fengyj@sem.tsinghua.edu.cn}.
#'
#' Filippo Palomba, Princeton University (maintainer). \email{fpalomba@princeton.edu}.
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}.
#'
#' @references
#' \itemize{
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}. Using synthetic controls: Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., and Titiunik, R.
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' scpi: Uncertainty Quantification for Synthetic Control Methods, \emph{arXiv}:2202.05984.}
#' \item{\href{https://arxiv.org/abs/2210.05026}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption, \emph{arXiv}:2210.05026.}
#' }
#'
#' @seealso \code{\link{scdata}}, \code{\link{scdataMulti}}, \code{\link{scest}}, \code{\link{scplot}}, \code{\link{scplotMulti}}
#'
#' @examples
#'
#' data <- scpi_germany
#'
#' df <- scdata(
#'   df = data, id.var = "country", time.var = "year",
#'   outcome.var = "gdp", period.pre = (1960:1990),
#'   period.post = (1991:2003), unit.tr = "West Germany",
#'   unit.co = setdiff(unique(data$country), "West Germany"),
#'   constant = TRUE, cointegrated.data = TRUE
#' )
#'
#' result <- scpi(df, w.constr = list(name = "simplex", Q = 1), cores = 1, sims = 10)
#' result <- scpi(df,
#'   w.constr = list(lb = 0, dir = "==", p = "L1", Q = 1),
#'   cores = 1, sims = 10
#' )
#'
#' @export

SCMs <- function(data,
                 outcome_models = c("None", "OLS", "Ridge", "Lasso", "AugSynth"),
                 feature_weights = c("uniform", "optimized"),
                 unit_weights = c("uniform", "optimized"),
                 w.constr = NULL,
                 V = "separate",
                 V.mat = NULL,
                 inference_type = "scpi",
                 solver = "ECOS",
                 P = NULL,
                 u.missp = TRUE,
                 u.sigma = "HC1",
                 u.order = 1,
                 u.lags = 0,
                 u.design = NULL,
                 u.alpha = 0.05,
                 e.method = "all",
                 e.order = 1,
                 e.lags = 0,
                 e.design = NULL,
                 e.alpha = 0.05,
                 sims = 200,
                 rho = NULL,
                 rho.max = 0.2,
                 lgapp = "generalized",
                 cores = 1,
                 plot = FALSE,
                 plot.name = NULL,
                 w.bounds = NULL,
                 e.bounds = NULL,
                 save.data = NULL,
                 verbose = TRUE) {
  if ((methods::is(data, "scdata") || methods::is(data, "scdataMulti")) == FALSE) {
    stop("data should be the object returned by running scdata or scdata_multi!")
  }

  if (methods::is(data, "scdata") == TRUE) {
    class.type <- "scpi_data"
  } else if (methods::is(data, "scdataMulti") == TRUE) {
    class.type <- "scpi_data_multi"
  }

  #############################################################################
  #############################################################################
  ## Estimation of synthetic weights
  # if (verbose) {
  #   cat("---------------------------------------------------------------\n")
  #   cat("Estimating Weights...\n")
  # }

  scm_model <- scest(data = data, w.constr = w.constr, V = V, V.mat = V.mat, solver = solver)

  sc.pred <- run_outcome_models(
    scm_model = scm_model,
    scm_data = data,
    treated_unit = data$specs$treated.units,
    outcome_models = outcome_models,
    period_post = data$specs$period.post,
    Y = data$Y.donors.post,
    Z0 = data$Y.donors,
    Z1 = data$Y.pre
  )

  #############################################################################
  #############################################################################
  ## Retrieve processed data from scest

  A <- sc.pred$data$A # Features of treated unit
  B <- sc.pred$data$B # Features of control units
  Z <- sc.pred$data$Z # B and C column-bind
  Y.donors <- data$Y.donors # Outcome variable of control units
  Y.donors.post <- data$Y.donors.post # Outcome variable of control units
  K <- sc.pred$data$specs$K # Number of covs for adjustment per feature
  KM <- sc.pred$data$specs$KM # Dimension of r (total number of covs for adj)
  J <- sc.pred$data$specs$J # Number of donors
  M <- sc.pred$data$specs$M # Number of features
  T0 <- sc.pred$data$specs$T0.features # Time periods used per feature
  T1 <- sc.pred$data$specs$T1.outcome # Number of out-of-sample periods
  features <- sc.pred$data$specs$features # Name of features
  constant <- sc.pred$data$specs$constant # Logical indicating whether a constant is included
  out.feat <- sc.pred$data$specs$out.in.features # Logical indicating whether the outcome variable is among features
  coig.data <- sc.pred$data$specs$cointegrated.data # Logical indicating whether B is cointegrated
  w.constr <- sc.pred$est.results$w.constr # Constraints on w
  V <- sc.pred$est.results$V # Weighting matrix
  w <- sc.pred$est.results$w # Estimated vector of weights
  r <- sc.pred$est.results$r # Estimated coefficients of covariates
  b <- sc.pred$est.results$b # w and r column-bind
  Y.post.fit <- sc.pred$est.results$Y.post.fit # Estimated post-treatment outcome for SC unit
  res <- sc.pred$est.results$res # Residuals from estimation
  outcome.var <- sc.pred$data$specs$outcome.var # name of outcome variable
  sc.effect <- sc.pred$data$specs$effect # Causal quantity of interest
  sparse.mat <- sc.pred$data$specs$sparse.matrices # Whether sparse matrices are involved or not
  om.model.preds <- sc.pred$yhats

  if (class.type == "scpi_data") {
    Jtot <- J
    KMI <- KM
    I <- 1
    T0.tot <- sum(T0) # Total number of observations used in estimation
    T0.M <- T0.tot
    T1.tot <- T1
    features <- list(features)
    out.feat <- list(out.feat)
    T0 <- list(T0)
    names(T0) <- sc.pred$data$specs$treated.units
  } else if (class.type == "scpi_data_multi") {
    J <- unlist(J)
    Jtot <- sum(J)
    KMI <- data$specs$KMI # total number of covariates used for adjustment
    I <- data$specs$I # number of treated units
    T0.M <- unlist(lapply(data$specs$T0.features, sum)) # observations per treated unit
    T0.tot <- sum(T0.M) # Total number of observations used in estimation
    T1.tot <- sum(unlist(T1)) # Total number of observations post-treatment
  }


  if (sc.effect == "unit") {
    T1.tot <- I
    T1 <- lapply(T1, function(x) 1)
  }

  # Check on P
  if (is.null(P) == TRUE) {
    P <- sc.pred$data$P # Matrix for out-of-sample prediction
  } else {
    # User-provided prediction matrix P (should be T1 by (J+KM))
    if (is.matrix(P) == FALSE) {
      stop("The object P should be a matrix!")
    }
    if (class.type == "scpi_data") {
      Prows <- T1
      Pcols <- J + KM
    } else if (class.type == "scpi_data_multi") {
      Prows <- T1.tot
      Pcols <- Jtot + KMI
    }

    if (nrow(P) != Prows) {
      stop(paste("The matrix P currently has", nrow(P), "rows when instead", Prows, "were expected
                 (i.e. the number of post-intervention periods)!"))
    }
    if (ncol(P) != Pcols) {
      stop(paste("The matrix P currently has", ncol(P), "columns when instead", Pcols, "were expected
                 (i.e. the size of the donor pool plus the number of covariates used in adjustment in the outcome equation)!"))
    }
  }

  if (!(e.method %in% c("gaussian", "ls", "qreg", "all"))) {
    stop("The object e.method should be one of 'gaussian', 'ls', 'qreg', or 'all'.")
  }

  if (!(lgapp %in% c("linear", "generalized"))) {
    stop("The object lgapp should be one of 'linear' or 'generalized'.")
  }

  # Check on u_sigma and user-provided bounds
  if (is.character(u.sigma) == FALSE) {
    stop("The object u.sigma should be of type character!!")
  } else {
    if (!(u.sigma %in% c("HC0", "HC1", "HC2", "HC3", "HC4"))) {
      stop("Supported variance estimators are 'HC0','HC1','HC2','HC3','HC4'.")
    }
  }

  if (is.null(w.bounds) == FALSE) {
    if (is.matrix(w.bounds) == FALSE) {
      stop("The object w.bounds should be a matrix!")
    }

    if (ncol(w.bounds) != 2) {
      stop("w.bounds should be a matrix with two columns: the first column for the
           lower bound, the second for the upper. In case you don't want to specify
           the lower or the upper bound just fill the specific column with NAs.")
    }

    if (nrow(w.bounds) != length(Y.post.fit)) {
      stop(paste0(
        "w.bounds should be a matrix with ", length(Y.post.fit),
        " rows (i.e. the number of post-intervention periods)."
      ))
    }
  }


  if (sims < 10) {
    stop("The number of simulations needs to be larger or equal than 10!")
  }

  if (is.null(e.bounds) == FALSE) {
    if (is.matrix(e.bounds) == FALSE) {
      stop("The object e.bounds should be a matrix!")
    }

    if (ncol(e.bounds) != 2) {
      stop("e.bounds should be a matrix with two columns: the first column for the
           lower bound, the second for the upper. In case you don't want to specify
           the lower or the upper bound just fill the specific column with NAs.")
    }

    if (nrow(e.bounds) != length(Y.post.fit)) {
      stop(paste(
        "e.bounds should be a matrix with ", length(Y.post.fit),
        " rows (i.e. the number of post-intervention periods)."
      ))
    }
  }

  # Check rho
  if (is.null(rho) == FALSE) {
    if (is.character(rho) == TRUE) {
      if (!(rho %in% c("type-1", "type-2", "type-3"))) {
        stop("When not a scalar, 'rho' must be 'type-1', 'type-2', or 'type-3'.")
      }
    }
  } else {
    rho <- "type-1"
  }

  # Check on number of cores
  if (is.null(cores) == FALSE) {
    n.cores <- parallel::detectCores(logical = TRUE)
    if (cores > n.cores) {
      stop(paste("You selected", cores, "cores, but only", n.cores, " cores are available on your machine!"))
    }
  } else {
    cores <- parallel::detectCores(logical = TRUE) - 1
    warning(paste("scpi is using", cores, "cores for estimation! You can adjust the number of cores with the 'cores' option."), immediate. = TRUE, call. = FALSE)
  }

  if (verbose) {
    if (class.type == "scpi_data") {
      constr.type <- w.constr[["name"]]
    } else if (class.type == "scpi_data_multi") {
      constr.type <- w.constr[[1]]$name
    }

    cat("Quantifying Uncertainty\n")

    # executionTime(T0.tot, Jtot, I, T1.tot, sims, cores, constr.type)
  }

  if (inference_type == "scpi") {
    inference.results <- inference_scpi(
      A = A,
      B = B,
      P = P,
      Z = Z,
      I = I,
      V = V,
      w = w,
      res = res,
      Y.donors = Y.donors,
      Y.post.fit = Y.post.fit,
      sc.pred = sc.pred,
      sc.effect = sc.effect,
      sparse.mat = sparse.mat,
      class.type = class.type,
      w.constr = w.constr,
      rho = rho,
      rho.max = rho.max,
      coig.data = coig.data,
      T0.M = T0.M,
      J = J,
      KM = KM,
      KMI = KMI,
      lgapp = lgapp,
      b = b,
      features = features,
      u.order = u.order,
      u.lags = u.lags,
      u.design = u.design,
      e.order = e.order,
      e.lags = e.lags,
      e.design = e.design,
      out.feat = out.feat,
      outcome.var = outcome.var,
      T0 = T0,
      T1 = T1,
      T1.tot = T1.tot,
      constant = constant,
      u.alpha = u.alpha,
      u.missp = u.missp,
      u.sigma = u.sigma,
      e.alpha = e.alpha,
      e.method = e.method,
      w.bounds = w.bounds,
      e.bounds = e.bounds,
      Jtot = Jtot,
      sims = sims,
      cores = cores,
      verbose = verbose
    )
  } else{
    
  }


  result <- list(
    data = sc.pred$data,
    est.results = sc.pred$est.results,
    inference.results = inference.results
  )

  class(result) <- "scpi"
  if (class.type == "scpi_data") {
    result$data$specs$class.type <- "scpi_scpi"
  } else if (class.type == "scpi_data_multi") {
    result$data$specs$class.type <- "scpi_scpi_multi"
  }

  #############################################################################
  #############################################################################
  ## Plot

  if (plot == TRUE) {
    if (is.null(plot.name) == FALSE) {
      fig.name <- plot.name
    } else {
      fig.name <- "scpi_default_plot"
    }

    if (class.type == "scpi_data") {
      scplot(
        result = result, fig.path = getwd(),
        fig.name = fig.name, fig.format = "png", save.data = save.data
      )
    } else {
      scplotMulti(result)
    }
  }

  return(result)
}
