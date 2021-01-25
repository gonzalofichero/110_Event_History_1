
# In the lecture video the data on female fruitflies were introduced and the Weibull distribution
# looked like a promising candidate model. Use these data (again remove the first two days and
# analyze only life spans after day 2) to do the following:


# Q1. Estimate the parameters a and b with maximum likelihood (numerical optimization).

flies <- read.table("flies.txt", header=T)

hist(flies$days,breaks=20)


# ------ remove 'infant shock' -----------
age2 <- flies$days[flies$days > 2] - 2   ## age2 is age in days after day 2
length(age2)   # now n=965

hist(age2,breaks=20)


# Loglikelihood function 
logL.W <- function(para, lifespans){
  a <- para[1]
  1
  b <- para[2]
  logL <- sum(log(dweibull(lifespans, shape=a, scale=b)))
  return(logL)
}

# Starting points
a.start <- 1
b.start <- 10
# check whether it works for starting values
logL.W(c(a.start, b.start), lifespans=age2)

# now optimize
mle.w <- optim(c(a.start, b.start), logL.W, lifespans=age2, control=list(fnscale=-1), 
               hessian=TRUE, method="BFGS")

# Plotting real data vs optimized model
hist(age2, freq=FALSE)
lines(0:50, dweibull(0:50, shape=mle.w$par[1], scale=mle.w$par[2]), col="red", lwd=2)


