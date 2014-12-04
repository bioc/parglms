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
 TRUE
}
