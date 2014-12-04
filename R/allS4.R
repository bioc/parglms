setOldClass("Registry")
setGeneric("parGLM", function(formula, store, ..., jobids) standardGeneric("parGLM"))
setMethod("parGLM", c("formula", "Registry"), function(formula, store, ..., jobids) {
 ndots = names(list(...))
 stopifnot(all(c("family", "binit", "tol", "maxit") %in% ndots))
 stopifnot(inherits(list(...)$family, "function"))
 stopifnot(inherits(list(...)$binit, "numeric"))
 stopifnot(inherits(list(...)$maxit, "numeric"))
 stopifnot(inherits(list(...)$tol, "numeric"))
 # try to check compatibility of length(binit) with features of formula?
 .parglm(formula, store, ..., jobids)
})
  
