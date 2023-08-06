#' Plot a time resolved mass spectrometry signal
#'
#' Plots the raw signal of a given isotope against time.
#' 
#' @param x an object of class \code{\link{timeresolved}} or
#'     \code{\link{PHdata}}
#' @param mass a string indicating the isotope of interest
#' @param label a string with the name of the run
#' @param run the run number
#' @param hop identifier for mass channel in Noblesse instruments (if
#'     \code{x} has class \code{WiscAr}).
#' @param tmin time stamp (in seconds) of 'time zero'
#' @param ... optional parameters
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' mMC <- loaddata(samplefile)
#' plot(mMC,"MD2-1a")
#' @rdname plot
#' @export
plot.timeresolved <- function(x,label=NULL,mass=NA,
                              run=1,hop='101',tmin=0,...){
    if (!is.null(label))
        run <- which(x$labels==label)-1
    if (length(run)!=1){
        print('invalid input into plot function')
        return(NA)
    }
    therun <- subset(x,i=run)
    nm <- length(x$masses)
    graphics::par(mfrow=rep(ceiling(sqrt(nm)),2))
    for (m in 1:nm){
        f <- fit(therun,tmin=tmin,mass=x$masses[m],returnfit=TRUE)
        p <- stats::predict.lm(f$fit,newdata=data.frame(tt=therun$thetime))
        graphics::plot(therun$thetime,therun$d[,m],type='p',
                       xlab='time',ylab=x$masses[m])
        graphics::lines(therun$thetime,p)
        graphics::points(0,f$intercepts,pch=22,bg='black')
    }
    f <- fit(therun,tmin=tmin,returnfit=TRUE)
    invisible(f$fit)
}
#' @rdname plot
#' @export
plot.WiscAr <- function(x,hop='101',...){
    plot(x[[hop]],...)
}
#' @examples
#' mPH <- loaddata(samplefile)
#' plot(mPH,"MD2-1a","Ar40")
#' @rdname plot
#' @export
plot.PHdata <- function(x,...){
    plot.timeresolved(x,...)
}

#' Plot a matrix with correlation coefficients
#'
#' Converts the covariance matrix to a correlation matrix and plots
#' this is a coloured image for visual inspection.
#' 
#' @param X a data structure (list) containing an item called `covmat' (covariance matrix)
#' @examples
#' graphics.off()
#' data(Melbourne)
#' plotcorr(Melbourne$X)
#' @export
plotcorr <- function(X){
    image.with.legend(z=stats::cov2cor(X$covmat),
                      color.palette=grDevices::heat.colors)
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
