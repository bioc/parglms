checkglms = function() {
 #
 # Dec 4 -- so far, linear/gaussian model demonstrated
 #
 require(MASS)
 data(anorexia)
 myr = makeRegistry("abc", file.dir=tempfile())
 chs = chunk(1:nrow(anorexia), n.chunks=18) # 4 recs/chunk
 f = function(x) anorexia[x,]
 batchMap(myr, f, chs)
 submitJobs(myr) # now getResult(myr,1) gives back a data.frame
 waitForJobs(myr) # simple dispersal
 pp = parGLM( Postwt ~ Treat + Prewt, myr,
   family=gaussian, binit = c(0,0,0,0), maxit=10, tol=.001 )
 theLM <- lm(Postwt~Treat+Prewt, data=anorexia)
 checkTrue(max(abs(pp$coefficients - coef(theLM))) < 1e-5 )
 require(sandwich)
 hc0 <- vcovHC(theLM, type="HC0")
 checkTrue(max(abs(pp$robust.variance - hc0)) < 1e-5)
#
# logistic regression
#
 library(datasets)
 data(infert)
 library(BatchJobs)
 library(parglms)
 myr = makeRegistry("abc", file.dir=tempfile())
 chs = chunk(1:nrow(infert), n.chunks=31) 
 f = function(x) infert[x,]
 batchMap(myr, f, chs)
 submitJobs(myr) 
 waitForJobs(myr) 
 m1 = glm(case~education+age+parity, data=infert,
      family=binomial)
 pp = parGLM(case~education+age+parity, myr,
    family=binomial, binit=rep(0, length(m1$coef)), tol=.001,
     maxit=50)
 checkTrue(max(abs(pp$coef-m1$coef)) < 1e-5)
 se1 = summary(m1)$coef[,2]
 sepp = sqrt(diag(pp$eff))
 checkTrue(max(abs(se1-sepp)) < 1e-5)

 closeAllConnections()
 TRUE
}
