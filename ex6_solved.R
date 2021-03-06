library(tidyverse)
set.seed(42)


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
plot(hazard)


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
         g.a = ((lifespan/b_hat)^a_hat) * (1/lifespan + (a_hat/lifespan) * log(lifespan/b_hat) ),
         g.b = (-1) * ((a_hat * a_hat) / (lifespan * b_hat)) * ((lifespan/b_hat)^a_hat),
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
lines(lifespan_ci$ci_low, col="red", lty = 3, lwd = 3)
lines(lifespan_ci$ci_high, col="red", lty = 3, lwd =3)



###################################################################################################
# Q6. Similarly, the estimated Weibull survival function S(t) has a standard
# error. Is this s.e. uniform or does it change with t? Can you explain the result?


# Hazard function
lifespan_surv <- as.data.frame(seq(1,60, by = 1))
names(lifespan_surv) <- "lifespan"

v11 <- V[1,1]
v21 <- V[2,1]
v12 <- V[1,2]
v22 <- V[2,2]

lifespan_surv %>% 
  mutate(survival = (pweibull(lifespan, shape = a_hat, scale = b_hat, lower.tail = FALSE)),
         a_hat = a_hat,
         b_hat = b_hat,
         g.a.s = (-1)*exp((-1)*((lifespan/b_hat)^a_hat)) * ((lifespan/b_hat)^a_hat) * log(lifespan/b_hat),
         g.b.s = (a_hat/b_hat)*exp((-1)*((lifespan/b_hat)^a_hat)) * ((lifespan/b_hat)^a_hat),
         v11 = v11,
         v21 = v21,
         v12 = v12,
         v22 = v22,
         se.q.s = sqrt(g.a.s * g.a.s * v11 + g.a.s * g.b.s * v21 + g.a.s * g.b.s * v12 + g.b.s * g.b.s * v22),
         ci_low = survival - 1.96 * se.q.s,
         ci_high = survival + 1.96 * se.q.s
  ) -> lifespan_surv_ci


# Plotting the final result of the delta method for h(t)
plot(lifespan_surv_ci$lifespan, lifespan_surv_ci$survival)
lines(lifespan_surv_ci$ci_low, col="red", lty = 1, lwd = 1)
lines(lifespan_surv_ci$ci_high, col="red", lty = 1, lwd = 1)


#############################################
# Alternative solution (R base) for Q6


# Create a function to do the matrix multiplication and just returning and atomic vector
# Loop for the desire vector of "time" and save the calculation inside a new vector to be used
# by R base


gradient_survival_weibull <- function(t){
  
  # Optimized parameters
  a_hat <-  mle.w$par[1]
  b_hat <-  mle.w$par[2]
  
  # Var-Cov matrix from optimization
  V <- V
  
  # Gradient for parameter a for S(t) for Weibull:
  g.a <- (-1)*exp((-1)*((t/b_hat)^a_hat)) * ((t/b_hat)^a_hat) * log(t/b_hat)
  # Gradient for parameter b for S(t) for Weibull:
  g.b <- (a_hat/b_hat)*exp((-1)*((t/b_hat)^a_hat)) * ((t/b_hat)^a_hat)
  
  # Putting it together
  gradient <- c(g.a, g.b)
  
  # Matrix multiplication:
  se.q <- sqrt( t(gradient) %*% V %*% gradient)
  
  # Getting back what we need
  return(se.q)
  
}

# Create empty vector to fill
se.q.gradient <- rep(NA,60)

# Loop through the possible values of t and apply the previous function to calculate
# standard errors for S(t) Weibull based on optimization
for (t in 0:60) {
  
  se.q.gradient[t] <- gradient_survival_weibull(t)
  
}


# Check
cbind(lifespan_surv_ci$lifespan, lifespan_surv_ci$se.q.s, se.q.gradient)
# All looks good :D