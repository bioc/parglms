### =========================================================================
### parglm tools
# idea is to allow parallel evaluation of contributions to 
# sufficient statistics
#
# use case is data sharded into a BatchJobs folder structure
# governed by a registry in a gQTLBase::ciseStore instance
#
# is it reasonable?  certainly if a BatchJobs result set exists
# and is to be analyzed, yes.
#
# what if we just have a very large data frame?  using BatchJobs
# to store it in shards allows persistent small-footprint retrievals
#
### -------------------------------------------------------------------------
###

job2df = function(store, i) {
 if (!is.null(store$extractor)) {
  stopifnot( all(names(formals(store$extractor))==c("store", "i")) )
  return(store$extractor(store, i))
  }
 loadResult(store, i)
}

getMu <- function (formula, store, i, beta, family) 
    as.numeric(family()$linkinv(getEta(formula, store, i, beta)))

getY <- function (formula, store, i) model.response(model.frame(formula, data=job2df(store, i)))

getX <- function (formula, store, i) 
{
    model.matrix(formula, job2df(store, i))
}

ri <- function(formula, store, i, beta, family) getY(formula, store,i) - 
     getMu(formula, store, i, beta, family)

### possibly move to C
getEta <- function (formula, store, i, beta) getX(formula, store, i) %*% beta

Di <- function (formula, store, i, beta, family) 
    as.numeric(family()$mu.eta(getEta(formula, store, i, beta))) * getX(formula, store, i)

#
# here we are going to return a vector
#
Vinv.i <- function(formula, store, i, beta, family) 
    1/family()$variance(getMu(formula, store, i, beta, family))

Gcomps <- function(formula, i, store, beta, family) {
    DD <- Di(formula, store, i, beta, family)
#    tDD <- t(DD)
    r_i <- ri(formula, store, i, beta, family)
    Vinv <- Vinv.i(formula, store, i, beta, family)
    val <- Vinv * DD
    val1 <- t(val) %*% DD
    val2 <- t(val) %*% r_i
    middle = crossprod( val * abs(r_i) ) # Xt diag(u_i^2) X so HC0 in Zeileis sandwich
    list(DtVDi=val1, DtVri=val2, DtVririVD=middle, rss=sum(r_i^2), totalN=length(r_i))
}

combi <- function(x) {
    xx <- x[[1]]
    for (i in 2:length(x))
        xx <- Map("+", xx, x[[i]]) 
    xx
}

.parglm <- function(formula, store, family, binit, maxit, tol, jobids, theCall) {
#
# idea is that Gcomps(formula, i, store, ...) computes elements on chunk i (no cluster structure)
#
    converged = FALSE
    beta <- binit
    del <- Inf
    curit <- 0
    robvar <- NA
    solve_DtVDi <- NA
    s2 = NA
    N = NA
    delcomp = NULL
#
# check compatibility of binit and X
#
    if (missing(jobids) & inherits(store, "Registry")) jobids = findDone(store)
    x1 = getX(formula, store, jobids[1])
    if (!(length(beta) == ncol(x1))) stop("length(binit) not compatible with X")
    while ( max(abs(del/beta)) > tol ) {
        if (maxit == 0) break
        if (curit > maxit) {
            message(paste0("NOTE: ", paste("maxit [", maxit, "] iterations exceeded")))
            break  # converged will be false
            }
#        res <- bplapply(jobids, function(ind) Gcomps(formula=formula, i=ind,
#		store=store, beta=beta, family=family))
        res <- foreach(ind=jobids) %dopar% { Gcomps(formula=formula, i=ind,
		store=store, beta=beta, family=family) }
        delcomp <- combi(res) 
        solve_DtVDi <- solve(delcomp$DtVDi)  # move to QR?
        del <- solve_DtVDi %*% delcomp$DtVri
        beta <- beta + del
        robvar <- solve_DtVDi %*% (delcomp$DtVririVD %*% solve_DtVDi) 
        if (options()$verbose) {
            print(paste("iter ", curit))
            print("beta:")
            print(beta)
            }
        curit <- curit + 1
    }
    N = delcomp$totalN
    s2 = delcomp$rss/(N - nrow(delcomp$DtVDi)) # sum(resids^2)/(N-nrow(delcomp[[1]]))
    if (maxit == 0) converged = NA
    else if (curit <= maxit) converged = TRUE
    fac = s2
    if (!(deparse(substitute(family))=="gaussian")) fac=1.0
    ans = list(coefficients=beta, eff.variance=fac*solve_DtVDi, robust.variance=robvar, s2=s2,
       niter = curit, converged=converged, formula=formula, N=N, theCall=theCall )
    class(ans) = c("parglm", "list")
    ans
}

summary.parglm = function(object, ...) {
#
# eventually we'll get an object
#
  co = object$coeff
  se = sqrt(diag(object$eff))
  robse = sqrt(diag(object$rob))
  z = object$coeff/se
  robz = object$coeff/robse
  data.frame(beta=co, s.e.=se, eff.z=z, rob.s.e.=robse, rob.z=robz)
}

predict.parglm = function (object, ...) { # (x, store, jobids) 
#
# objective is to apply coefficients built from one
# set of jobids to data in a complementary set of jobids
#
    inextras = list(...)
    store = inextras$store
    if (is.null(store)) stop("please supply registry as value of argument 'store'")
    jobids = inextras$jobids
    if (is.null(jobids)) stop("please supply vector of jobids to supply data for prediction")
    f = store$extractor
    if (is.null(f)) f = loadResult
    #bplapply(jobids, function(z) {
    foreach(z = jobids) %dopar% {
        tmp = f(store, z)
        fr = model.frame(object$formula, tmp)
        obs = model.extract(fr, "response")
        xx = model.matrix(object$formula, data = tmp)
        xx %*% object$coefficients
    }   #)
}

print.parglm = function(x, ...) {
 cat("parGLM result.  the call was:\n")
 show(x$theCall)
 if (isTRUE(x$converged)) cat("converged in", x$niter, "iterations.\n")
 else cat("failed to converge in", x$niter, "iterations.\n")
 show(summary(x))
}
