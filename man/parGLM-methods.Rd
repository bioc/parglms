\name{parGLM-methods}
\docType{methods}
\alias{parGLM-methods}
\alias{parGLM,formula,Registry-method}
\alias{parGLM}
\title{ fit GLM-like models with parallelized contributions to sufficient statistics }
\description{
This package addresses the problem of fitting GLM-like models in a scalable way,
recognizing that data may be dispersed, with chunks processed in parallel,
to create low-dimensional summaries from which model fits may be constructed.
}
\section{Methods}{
\describe{

\item{\code{signature(formula = "formula", store = "Registry")}}{
The model data are assumed to lie in the \code{file.dir/jobs/*}
folders, with \code{file.dir} defined in the \code{store}, which is
an instance of \code{\link[BatchJobs]{Registry}}.

Additional arguments must be supplied:
\describe{
   \item{family}{a function that serves as
a family for stats::glm}
   \item{binit}{a vector of initial values for regression
parameter estimation, must conform to expectations of \code{formula}} 
  \item{maxit}{an integer giving the maximum number of iterations allowed}
  \item{tol}{a numeric giving the tolerance criterion}
}
Failure to specify these triggers a fatal error.

The Registry instance can be modified to include a list element
'extractor'.  This must be a function with arguments \code{store}, and
code{i}.  The standard extraction function is 

\code{function(store, i) loadResult(store, i)}

It must return a data frame, conformant with the expectations of \code{formula}.
Limited checking is performed.
}
}}
\examples{
if (require(MASS) & require(BatchJobs)) {
 data(anorexia)  # N = 72
 myr = makeRegistry("abc", file.dir=tempfile())
 chs = chunk(1:nrow(anorexia), n.chunks=18) # 4 recs/chunk
 f = function(x) anorexia[x,]
 batchMap(myr, f, chs)
 submitJobs(myr) # now getResult(myr,1) gives back a data.frame
 waitForJobs(myr) # simple dispersal
 pp = parGLM( Postwt ~ Treat + Prewt, myr,
   family=gaussian, binit = c(0,0,0,0), maxit=10, tol=.001 )
 print(summary(theLM <- lm(Postwt~Treat+Prewt, data=anorexia)))
 print(pp$coefficients - coef(theLM))
 if (require(sandwich)) {
   hc0 <- vcovHC(theLM, type="HC0")
   print(pp$robust.variance - hc0)
   }
 }
}
\keyword{methods}
\keyword{ modeling }
