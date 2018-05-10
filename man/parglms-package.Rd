\name{parglms-package}
\alias{parglms-package}
\alias{parglms}
\docType{package}
\title{
\packageTitle{parglms}
}
\description{
\packageDescription{parglms}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{parglms}
\packageIndices{parglms}

In version 0.0.0 we established an approach to fitting GLM from
data that have been persistently dispersed and managed by 
a \code{\link{Registry}}.
}
\author{
\packageAuthor{parglms}

Maintainer: \packageMaintainer{parglms}
}
\references{
This package shares an objective with the \code{\link{bigglm}}
methods of \code{biglm}.  In \code{bigglm}, a small-RAM-footprint algorithm
is employed, with sequential chunking to update statistics in each iteration.
In \code{parGLM} the footprint is likewise controllable, but statistics
in each iteration are evaluated in parallel over chunks.
}
\keyword{ package }
%\seealso{
%\code{\link[bigglm:-package]{<pkg>}} ~~
%}
\examples{
showMethods("parGLM")
}
