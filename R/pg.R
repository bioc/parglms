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

#job2dfClo = function(tx=force) function(store, i) tx( loadResult(store, i) )

#job2df = job2dfClo()
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
    list(DtVDi=val1, DtVri=val2, DtVririVD=middle, residi=r_i)
}

combi <- function(x) {
    kp = c("DtVDi", "DtVri", "DtVririVD") # leave residi alone
    xx <- x[[1]][kp]
    for (i in 2:length(x))
        xx <- Map("+", xx, x[[i]][kp]) 
    xx
}

.parglm <- function(formula, store, family, binit, maxit, tol, jobids) {
#
# idea is that Gcomps(formula, i, store, ...) computes elements on chunk i (no cluster structure)
#
    beta <- binit
    del <- Inf
    curit <- 0
    robvar <- NA
#
# check compatibility of binit and X
#
    if (missing(jobids) & inherits(store, "Registry")) jobids = findDone(store)
    x1 = getX(formula, store, jobids[1])
    if (!(length(beta) == ncol(x1))) stop("length(binit) not compatible with X")
    while (max(abs(del/beta)) > tol ) {
        res <- bplapply(jobids, function(ind) Gcomps(formula=formula, i=ind,
		store=store, beta=beta, family=family))
        delcomp <- combi(res) 
        resids = unlist(lapply(res, function(x) x$residi))  # vector will be long -- ff?
        N = length(resids)
        s2 = sum(resids^2)/(N-nrow(delcomp[[1]]))
        solve_DtVDi <- solve(delcomp[[1]])  # move to QR?
        del <- solve_DtVDi %*% delcomp[[2]]
        beta <- beta + del
        robvar <- solve_DtVDi %*% (delcomp[[3]] %*% solve_DtVDi) 
        if (options()$verbose) {
            print(paste("iter ", curit))
            print("beta:")
            print(beta)
            }
        curit <- curit + 1
        if (curit > maxit) 
            stop(paste("maxit [", maxit, "] iterations exceeded"))
    }
    list(coefficients=beta, eff.variance=s2*solve_DtVDi, robust.variance=robvar, s2=s2,
       niter = curit-1 )
      
}

