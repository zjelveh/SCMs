#' @title Plot Synthetic Control Point Estimates and Prediction Interval
#'
#' @description Plot treated and synthetic time series from a fitted
#' \code{\link{scest}} object.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are \link{scdata} for data preparation, \link{scest} for
#' point estimation, and \link{inference_sc} for inference procedures.
#'
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#'
#' @param result A class \code{"scest"} object from \code{\link{scest}}.
#' @param fig.path a string indicating the path where the plot(s) should be saved.
#' @param fig.name a string indicating the name of the plot(s). If multiple plots will be saved the command automatically
#' generates a numeric suffix to avoid overwriting them.
#' @param fig.format a string indicating the format in which the plot(s) should be saved.
#' @param col.treated a string specifying the color for the treated unit series. Find the full list at \url{http://sape.inf.usi.ch/quick-reference/ggplot2/colour}.
#' @param col.synth a string specifying the color for the synthetic unit series. Find the full list at \url{http://sape.inf.usi.ch/quick-reference/ggplot2/colour}. 
#' @param label.xy a character list with two elements indicating the name of the axes
#' (eg. label.xy = list(x.lab = "Year", y.lab = "GDP growth (%)")).
#' @param plot.range a numeric array indicating the time range of the plot(s).
#' @param x.ticks a numeric list containing the location of the ticks on the x axis.
#' @param event.label a list containing a character object ('lab') indicating the label of the event and
#' a numeric object indicating the height of the label in the plot.
#' @param plot.specs a list containing some specifics to be passed to ggsave (eg. img.width, img.height, dpi)
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used
#' to produce the plot.
#'
#' @return
#' \item{plots}{a list containing standard ggplot object(s) that can be used for further customization.}
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
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}. Using synthetic controls:
#' Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' }
#'
#' @seealso \code{\link{scdata}}, \code{\link{scest}}, \code{\link{inference_sc}}
#'
#' @examples
#' \dontrun{
#' # Create example data
#' df <- data.frame(
#'   unit = rep(c("treated", "control1", "control2"), each = 20),
#'   time = rep(1:20, 3),
#'   outcome = c(rnorm(20, mean = 1:20), rnorm(20, mean = 1:20), rnorm(20, mean = 1:20))
#' )
#'
#' # Prepare data for SCM
#' scm_data <- scdata(df = df, id.var = "unit", time.var = "time",
#'                    outcome.var = "outcome", period.pre = 1:10,
#'                    period.post = 11:20, unit.tr = "treated",
#'                    unit.co = c("control1", "control2"))
#'
#' # Estimate synthetic control
#' result <- scest(scm_data, w.constr = list(name = "simplex"))
#'
#' # Plot results
#' scplot(result)
#' }
#'
#' @export
#'
#'

scplot  <- function(result, fig.path = NULL, fig.name = NULL, fig.format = "png",
                    col.treated = "black", col.synth = "mediumblue", label.xy = NULL, plot.range = NULL, 
                    x.ticks = NULL, event.label = NULL, plot.specs = NULL, save.data = NULL) {

  if (!is.null(fig.path)) {
    if (is.character(fig.path) == FALSE) stop("The object 'fig.path' should be a character!")
  }

  if (!is.null(fig.name)) {
    if (is.character(fig.name) == FALSE) stop("The object 'fig.name' should be a character!")
  }

  if (is.character(fig.format) == FALSE) stop("The object 'fig.format' should be a character!")

  if (!is.null(save.data)) {
    if (is.character(save.data) == FALSE) stop("The object 'save.data' should be a character!")
  }

  if (is.null(label.xy) == FALSE){
    if (is.list(label.xy) == FALSE) {
      stop("label.xy should be a list of two elements! (eg. label.xy = list(x.lab = 'Year', y.lab = 'GDP')")
    }

    if (!all(names(label.xy) %in% c("x.lab", "y.lab"))) {
      stop("label.xy should be a list of two elements named 'x.lab' ad 'y.lab'!")
    }

    if (!all(sapply(c(label.xy$x.lab, label.xy$y.lab), is.character))) {
      stop("label.xy should be a list of two character elements!")
    }
  }

  if (is.null(label.xy) == TRUE) {
    x.lab <- "Time"
    y.lab <- "Outcome Variable"
  } else {
    x.lab <- label.xy$x.lab
    y.lab <- label.xy$y.lab
  }

  if (is.null(event.label) == FALSE){
    if (is.list(event.label) == FALSE) {
      stop("event.label should be a list of two elements! (eg. event.label = list(lab = 'Event', y.lab = 10)")
    }

    if (!all(names(event.label) %in% c("lab", "height"))) {
      stop("event.label should be a list of two elements named 'lab' ad 'height'!")
    }

    if (!all(is.character(event.label$lab), is.numeric(event.label$height))) {
      stop("event.label should be a list of two elements: one numeric ('height'), one character ('lab')!")
    }
  } 

  if (is.null(plot.specs) == FALSE){
    if (is.list(plot.specs) == FALSE) {
      stop("plot.specs should be a list of elements! (eg. plot.specs = list(img.height = 4.5, img.width = 6, dpi = 1000)")
    }

    if (!all(names(plot.specs) %in% c("img.height", "img.width", "dpi"))) {
      stop("plot.specs should be a list of elements named 'img.height', 'img.width', 'dpi'!")
    }

    img.width  <- plot.specs$img.width
    img.height <- plot.specs$img.height
    dpi        <- plot.specs$dpi
  } else {
    img.width  <- 6
    img.height <- 4.5
    dpi        <- 1000

  }

  # Some plot specifics used if the plot is saved
  if (any(is.null(fig.path), is.null(fig.name))) {
      save.plot  <- FALSE
  } else {
      save.plot  <- TRUE
  }

  if (!inherits(result, "scest")) {
    stop("The object 'result' should be a class 'scest' object produced by scest().")
  }

  if (save.plot == TRUE) {
    if (!fig.format %in% c("eps","ps","tex","pdf","jpeg","tiff","png","bmp","svg","wmf")) {
      stop("The specified format is not valid. See ?ggsave to check valid formats.")
    } else {
      plot.fitted <- paste(fig.path,"/",fig.name,".",fig.format, sep = "")
    }
  }

  y.fit   <- rbind(result$est.results$Y.pre.fit, result$est.results$Y.post.fit)
  y.act   <- rbind(result$data$Y.pre, result$data$Y.post)


  period.pre  <- result$data$specs$period.pre
  period.post <- result$data$specs$period.post
  T0 <- period.pre[length(period.pre)] # intercept

  if (is.null(plot.range)) {
    plot.range <- c(period.pre, period.post)
  }


  # Actual data
  dat    <- data.frame(t     = c(period.pre, period.post),
                       Y.act = c(y.act),
                       sname = "Treated")

  # Fill with NAs Y.fit and confidence bounds where missing
  Y.fit.na  <- matrix(NA, nrow = length(c(period.pre, period.post)))
  names <- strsplit(rownames(y.fit), "\\.")
  rnames <- unlist(lapply(names, "[[", 2))
  not.missing.plot <- as.character(c(period.pre,period.post)) %in% rnames

  Y.fit.na[not.missing.plot, 1]  <- y.fit

  # Synthetic unit data
  dat.sc <- data.frame(t        = c(period.pre, period.post),
                       Y.sc     = Y.fit.na,
                       sname    = "SC Unit")

  if (is.null(x.ticks)) {
    x.ticks <- c(seq(plot.range[1], plot.range[length(plot.range)], length.out = 5), T0)
    x.ticks <- round(unique(x.ticks))
  }

  if (is.null(event.label)) {
    event.lab <- ""
    event.lab.height <- dat.sc$Y.sc[dat.sc$t == T0]
  } else {
    event.lab <- paste("\n", event.label$lab, sep = "")
    event.lab.height <- event.label$height
  }


  dat.plot    <- dat[dat[,'t'] %in% plot.range, ]
  dat.sc.plot <- dat.sc[dat.sc[,'t'] %in% plot.range, ]

  plotdf <- dplyr::left_join(dat.plot,dat.sc.plot, by = 't')

  ## Plot specs
  plot <- ggplot() + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_text(vjust = 0.5, hjust=1)) +
          labs(x = x.lab, y = y.lab) +
          theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_blank(),
                legend.background = element_rect(fill = "white", color = "black"))

  ## Series to plot
  plot  <- plot + geom_line( data = plotdf, aes(x = t, y = .data$Y.act, colour = .data$sname.x), linetype = 'solid') +
                  geom_point(data = plotdf, aes(x = t, y = .data$Y.act, colour = .data$sname.x), shape = 1) +
                  geom_line( data = plotdf, aes(x = t, y = .data$Y.sc,  colour = .data$sname.y), linetype = 'dashed') +
                  geom_point(data = plotdf, aes(x = t, y = .data$Y.sc,  colour = .data$sname.y), shape = 19) +
                  geom_vline(xintercept = T0, linetype = "dashed") +
                  geom_text(aes(x = T0, label = event.lab, y = event.lab.height), angle = 90, size = 4) +
                  scale_x_continuous(breaks = x.ticks) +
                  scale_color_manual(name = "", values = c(col.synth, col.treated),
                                     labels = c("Synthetic Control", "Treated"),
                                     guide = guide_legend(override.aes = list(
                                       linetype = c('dashed','solid'), shape = c(19, 1)))) +
                  ggtitle("Synthetic Control Prediction")

  if (save.plot == TRUE) {
    suppressWarnings(ggsave(filename = plot.fitted, plot = plot, width = img.width, height = img.height, dpi = dpi))
    cat("Plot saved at '", fig.path,"'\n", sep = "")
    cat("File name: ", plot.fitted,"\n", sep = "")
  }


  ## Return plot object to be modified by the user if needed
  plots <- list('plot_out' = plot)

  ## Save data to reproduce plot
  if (!is.null(save.data)) {
    save(plotdf, file = paste0(save.data, ".RData"))
  }

  return(plots)
}
