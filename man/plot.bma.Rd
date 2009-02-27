\name{plot.bma}
\alias{plot.bma}
\title{Plot Diagnostics for an blm Object}
\description{ Four plots (selectable by 'which') are currently
  available: a plot
  of residuals against fitted values, Cumulative Model Probabilities,
  log marginal likelihoods versus model dimension, and marginal
  inclusion probabilities.
}
\usage{\method{plot}{bma}(x, which=c(1:4),caption = c("Residuals vs Fitted",
                 "Model Probabilities", "Model Complexity",
                 "Inclusion Probabilities"),
                 panel = if (add.smooth) panel.smooth
                 else points, sub.caption = NULL, main = "", ask =
                 prod(par("mfcol")) < length(which) &&
                 dev.interactive(), ..., id.n = 3, labels.id =
                 names(residuals(x)), cex.id = 0.75, add.smooth =
                 getOption("add.smooth"), label.pos = c(4, 2))
}
\arguments{
  \item{x}{\code{blm} object, typically result of 'blm'}
   \item{which}{ if a subset of the plots is required, specify a subset of the
          numbers '1:4'}.
 \item{caption}{captions to appear above the plots}
 \item{panel}{panel function.  The useful alternative to 'points',
          'panel.smooth' can be chosen by 'add.smooth = TRUE'}
\item{sub.caption}{ common title-above figures if there are multiple; used as
          'sub' (s.'title') otherwise.  If 'NULL', as by default, a
          possible shortened version of \code{deparse(x$call)} is used}
\item{main}{title to each plot-in addition to the above 'caption'}
\item{ask}{logical; if 'TRUE', the user is asked before each plot, see
          'par(ask=.)'}
\item{...}{ other parameters to be passed through to plotting functions}
\item{id.n}{number of points to be labelled in each plot, starting with
          the most extreme}
\item{labels.id}{ vector of labels, from which the labels for extreme points
          will be chosen.  'NULL' uses observation numbers}
\item{cex.id}{ magnification of point labels.}
\item{add.smooth}{ logical indicating if a smoother should be added to most
          plots; see also 'panel' above}
\item{label.pos}{positioning of labels, for the left half and right half of
          the graph respectively, for plots 1-3}
}
\details{
}
\author{Merlise Clyde, based on plot.lm by John Maindonald and Martin Maechler}
\references{
}
\seealso{ \code{\link{plot.coef.bma}} and \code{\link{image.bma}}.
}
\examples{

}

\keyword{regression}
