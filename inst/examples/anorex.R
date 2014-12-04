
library(BiocParallel)
library(MASS)
data(anorexia)
library(parglms)
library(BatchJobs)
myr = makeRegistry("abc", file.dir=tempfile())
chs = chunk(1:nrow(anorexia), n.chu=18)
f = function(x) anorexia[x,]
batchMap(myr, f, chs)
submitJobs(myr) # now getResult(myr,1) gives back a data.frame

pp = parGLM( Postwt ~ Treat + Prewt, myr, 
   family=gaussian, binit = c(0,0,0,0), maxit=10,
   sandwich=TRUE, tol=.001 )
print(pp)
print(summary(theLM <- lm(Postwt~Treat+Prewt, data=anorexia)))
library(sandwich)
print(hc0 <- vcovHC(theLM, type="HC0"))

