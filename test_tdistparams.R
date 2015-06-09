require (fBasics)
set.seed (3835285)
y <- 2 + 3 * rt (n=300, df=5.3)
result <- tFit (y, doplot=FALSE)
result <- attributes (result)
print (result$fit$estimate)


source ("./tdistparams.R")
# Let us try `get.tparams.only` on our simulated data set $y$:

my.estimate <- get.tparams.only (y, 2.1, 20)
str <- sprintf ("nu_hat = %.2f, c_hat = %.2f, s_hat = %.2f",
	my.estimate$nu, my.estimate$Location, my.estimate$Scale)
print (str)



ptm <- proc.time()
CL <- 0.95
result <- get.tparams.withCI (y, 2.1, 20, CL = CL) 
proc.time() - ptm
str <- sprintf ("nu = %.2f, %.0f%% CI = [%.2f, %.2f]",
	result$nu, (100*CL), result$nu.CI[1], result$nu.CI[2])
print (str)
str <- sprintf ("Location = %.2f, %.0f%% CI = [%.2f, %.2f]",
	result$Location, (100*CL), result$Location.CI[1], result$Location.CI[2])
print (str)
str <- sprintf ("Scale = %.2f, %.0f%% CI = [%.2f, %.2f]",
	result$Scale, (100*CL), result$Scale.CI[1], result$Scale.CI[2])
print (str)



# References

