\name{predict.bma}
\alias{predict.bma}
\alias{predict}
\title{Prediction Method for an object of class BMA}
\description{Predictions under model averaging from a BMA object
}
\usage{\method{predict}{bma}(object, newdata, top=NULL, ...)}

\arguments{
  \item{object}{An object of class BMA, created by \code{bas}}
  \item{newdata}{new data for predictions}
  \item{top}{Use only the top M models, based on posterior
    probabilities}
  \item{...}{optional extra arguments}
}
\value{a list of
  \item{Ybma}{predictions using BMA}
  \item{Ypred}{matrix of predictions under each model}
  \item{best}{index of top models included}
}
\details{
}
\references{
}
\seealso{\code{\link{bas}}, \code{\link{fitted.bma}}
}
\examples{
\dontrun{data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4, alpha=13, prior="g-prior", initprobs="Uniform")
predict(hald.gprior, hald.gprior$X, top=5)
}}
\author{Merlise Clyde}
\keyword{regression}
