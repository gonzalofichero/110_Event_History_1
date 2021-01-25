library(tidyverse)


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




# Optimized Weibull hazard
hazard <- dweibull(1:60, shape = a_hat, scale = b_hat)/(1-pweibull(1:60, shape = a_hat, scale = b_hat))
plot(log(hazard))


################################################################################
# Q3. Plot the resulting survival function from the estimated Weibull distribution and add the
# empirical survival function in the same figure. From the figure, do you think that the
# Weibull distribution is a good model?


# Plotting real data vs optimized model: f(t)
# And then pass to S(t)
hist(age2, freq=FALSE, ylim = c(0,0.06))
lines(density(age2, adjust =2), col="black", lwd = 2)
lines(0:60, dweibull(0:60, shape=mle.w$par[1], scale=mle.w$par[2]), col="red", lwd=2)



# Empirical Survival function

age.death <- sort(unique(age2))   # unique observed ages at death
no.age.death <- table(age2)       # counts


#  what is empirical survival function? Proportion surviving age t
n <- sum(no.age.death)
F.emp <- cumsum(no.age.death)/n   # empirical cdf 
S.emp <- 1 - F.emp                # S = 1 - F
plot(age.death, S.emp, type="s")  # plot as step function
lines(0:50, pweibull(0:50, shape=mle.w$par[1], scale=mle.w$par[2], lower.tail = FALSE), col="red", lwd=2) #add optimized Weibull
title(main="Survival function (after removing infant deaths)", font.main=1)



#################################################################################################
# Q4. What are the standard errors for a and b? Please give 95% confidence intervals for the two
# parameters. Are the two estimates correlated?

# Build variance-covariance from Hessian 
V <- solve(-mle.w$hessian)      # inverse of negative Hessian

#  square root of diagonal elements are the s.e. of a and b
sqrt(diag(V))
# parameter estimates correlated?
cov2cor(V)

# CI for a and b (assuming 1-alpha=95%)
se.a <- sqrt(diag(V))[1]
se.b <- sqrt(diag(V))[2]

# Confidence intervals for both a and b (estimated) for Weibull function
c(a_hat - 1.96*se.a , a_hat + 1.96*se.a) 
c(b_hat - 1.96*se.b , b_hat + 1.96*se.b) 


###################################################################################################
# Q5. As a and b have some uncertainty, also the estimated hazard, which depends on a and b,
# has uncertainty, that is, it has a standard error. Derive the standard error for h(t) and plot a
# 95% confidence interval around the hazard (see item 2).


# Hazard function
hazard

lifespan <- as.data.frame(seq(1,60, by = 1))
names(lifespan) <- "lifespan"

v11 <- V[1,1]
v21 <- V[2,1]
v12 <- V[1,2]
v22 <- V[2,2]

lifespan %>% 
  mutate(hazard = dweibull(lifespan, shape = a_hat, scale = b_hat)/(1-pweibull(lifespan, shape = a_hat, scale = b_hat)),
         a_hat = a_hat,
         b_hat = b_hat,
         g.a = (1/b_hat) * ((lifespan/b_hat)^(a_hat)) * (1 + (a_hat/b_hat) * (log(lifespan/b_hat))),
         g.b = ((-1)*(1/(b_hat^2))) * (a_hat * ((lifespan/b_hat)^(a_hat-1)) + (a_hat - 1 ) * ((lifespan/b_hat)^(a_hat-2)) * lifespan),
         v11 = v11,
         v21 = v21,
         v12 = v12,
         v22 = v22,
         se.q = sqrt(g.a * g.a * v11 + g.a * g.b * v21 + g.a * g.b * v12 + g.b * g.b * v22),
         ci_low = hazard - 1.96 * se.q,
         ci_high = hazard + 1.96 * se.q
         ) -> lifespan_ci
  
# Plotting the final result of the delta method for h(t)
plot(lifespan_ci$lifespan, lifespan_ci$hazard)
lines(lifespan_ci$ci_low)
lines(lifespan_ci$ci_high)

