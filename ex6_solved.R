
# In the lecture video the data on female fruitflies were introduced and the Weibull distribution
# looked like a promising candidate model. Use these data (again remove the first two days and
# analyze only life spans after day 2) to do the following:


##############################################################################################
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


# Optimized a is
a_hat <- mle.w$par[1]

# Optimized b is
b_hat <- mle.w$par[2]


################################################################################
# Q2. Plot the resulting hazard rate h(t) and interpret briefly

# Plotting real data vs optimized model
hist(age2, freq=FALSE, ylim = c(0,0.06))
lines(density(age2, adjust =2), col="black", lwd = 2)
lines(0:50, dweibull(0:50, shape=mle.w$par[1], scale=mle.w$par[2]), col="red", lwd=2)


################################################################################
# Q3. Plot the resulting survival function from the estimated Weibull distribution and add the
# empirical survival function in the same figure. From the figure, do you think that the
# Weibull distribution is a good model?


# Empirical Survival function

age.death <- sort(unique(age2))   # unique observed ages at death
no.age.death <- table(age2)       # counts

#  what is empirical survival function? Proportion surviving age t
n <- sum(no.age.death)
F.emp <- cumsum(no.age.death)/n   # empirical cdf 
S.emp <- 1 - F.emp                # S = 1 - F
plot(age.death, S.emp, type="s")  # plot as step function
lines(0:50, rev(pweibull(0:50, shape=mle.w$par[1], scale=mle.w$par[2])), col="red", lwd=2) #add optimized Weibull
title(main="Survival function (after removing infant deaths)", font.main=1)



################################################################################
# Q4. What are the standard errors for a and b? Please give 95% confidence intervals for the two
# parameters. Are the two estimates correlated?

# Build variance-covariance from Hessian 
V <- solve(-mle.w$hessian)      # inverse of negative Hessian

#  square root of diagonal elements are the s.e. of a and b
sqrt(diag(V))
# parameter estimates correlated?
cov2cor(V)

# CI for a and b (assuming 1-alpha=95%)
se.a <- sqrt(diag(V)) [1]
se.b <- sqrt(diag(V)) [2]

# Confidence intervals for both a and b (estimated) for Weibull function
c(a_hat - 1.96*se.a , a_hat + 1.96*se.a) 
c(b_hat - 1.96*se.b , b_hat + 1.96*se.b) 


# Q5. As a and b have some uncertainty, also the estimated hazard, which depends on a and b,
# has uncertainty, that is, it has a standard error. Derive the standard error for h(t) and plot a
# 95% confidence interval around the hazard (see item 2.).


# --- estimate median life span 
q.hat <- 1/b.hat * log(1 - b.hat*log(0.5)/a.hat )
q.hat  # Median age is q.hat+50

plot(0:60, 1-exp(- a.hat/b.hat *(exp(b.hat*(0:60))-1)), type="l", lwd=2,
     xlab="x", ylab="F(x)")
abline(h=0.5, lty=2)


C <- b.hat * log(0.5) / a.hat
# now gradient
g.a <- 1/(a.hat * b.hat) * C/(1-C)
g.b <- -1/b.hat^2 *(log(1-C) + C/(1-C))
gradient <- c(g.a, g.b)

se.q <- sqrt( t(gradient) %*% V %*% gradient)
se.q      # is 1x1 

#    confidence interval for q 
c(q.hat - 1.96 * se.q , q.hat + 1.96 * se.q)
# let's add 50 to obtain original ages
c(q.hat - 1.96 * se.q , q.hat + 1.96 * se.q)  + 50
