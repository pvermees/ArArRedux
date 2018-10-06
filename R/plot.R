#' Plot a time resolved mass spectrometry signal
#'
#' Plots the raw signal of a given isotope against time.
#' 
#' @param x an object of class \code{\link{timeresolved}} or
#' \code{\link{PHdata}}
#' @param mass a string indicating the isotope of interest
#' @param label a string with the name of the run
#' @param run the run number
#' @param ... optional parameters
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' mMC <- loaddata(samplefile,masses)
#' plot(mMC,"MD2-1a","Ar40")
#' @rdname plot
#' @export
plot.timeresolved <- function(x,mass,label=NULL,run=1,...){
    timeresolvedplot(x,mass,label=label,run=run,...)
}
#' @examples
#' mPH <- loaddata(samplefile,masses,PH=TRUE)
#' plot(mPH,"MD2-1a","Ar40")
#' @rdname plot
#' @export
plot.PHdata <- function(x,mass,label=NULL,run=1,...){
    plot.timeresolved(x$signals[[mass]],mass,label=label,run=run,...)
}
timeresolvedplot <- function(x,...){ UseMethod("timeresolvedplot",x) }
timeresolvedplot.default <- function(x,...){stop()}

#' Plot a matrix with correlation coefficients
#'
#' Converts the covariance matrix to a correlation matrix and plots
#' this is a coloured image for visual inspection.
#' 
#' @param X a data structure (list) containing an item called `covmat' (covariance matrix)
#' @examples
#' data(Melbourne)
#' plotcorr(Melbourne$X)
#' @export
plotcorr <- function(X){
    image.with.legend(z=stats::cov2cor(X$covmat),color.palette=grDevices::heat.colors)
}

# modified version of filled.contour with ".filled.contour" part replaced with "image"
# function. Note that the color palette is a flipped heat.colors rather than cm.colors
image.with.legend <- function (x = seq(1, nrow(z), length.out = nrow(z)), y = seq(1, 
    ncol(z), length.out=nrow(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = grDevices::heat.colors, 
    col = rev(color.palette(length(levels) - 1)), plot.title, plot.axes,
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) {
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(1, nrow(z), length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- graphics::par(c("mar", "las", "mfrow")))$mar
    on.exit(graphics::par(par.orig))
    w <- (3 + mar.orig[2L]) * graphics::par("csi") * 2.54
    graphics::layout(matrix(c(2, 1), ncol = 2L), widths = c(1, graphics::lcm(w)))
    graphics::par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    graphics::par(mar = mar)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = range(levels),
                xaxs = "i", yaxs = "i")
    graphics::rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
        if (axes) 
            graphics::axis(4)
    }
    else key.axes
    graphics::box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    graphics::par(mar = mar)
    graphics::image(x,y,z,col=col,xlab="",ylab="")
    if (missing(plot.axes)) {
        if (axes) {
            graphics::title(main = "", xlab = "", ylab = "")
            graphics::Axis(x, side = 1)
            graphics::Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        graphics::box()
    if (missing(plot.title)) 
        graphics::title(...)
    else plot.title
    invisible()
}
