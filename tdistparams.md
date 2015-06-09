---
title: "Estimating the parameters of a t-distributed random variable"
author: "Dr. Hossain Pezeshki"
date: '2015-06-06'
output: html_document
css: style.css
bibliography: bibliography.bib
csl: ieee.csl
---

# Prompt #
One notices that the `tFit{fBasics}` function @fBasics in R version 3.2.0 @theRlang
does not work quite right.
In fact, let us generate one hundred independent samples of a
$t_{\nu=5.3}$ random variable with location
parameter $2$ and scale parameter $3$, and apply `tFit` to the resulting data as follows:

```r
require (fBasics)
set.seed (3835285)
y <- 2 + 3 * rt (n=300, df=5.3)
result <- tFit (y, doplot=FALSE)
result <- attributes (result)
print (result$fit$estimate)
```

```
##             df 
## 0.829248854222
```
Namely the degrees-of-freedom is estimated as $0.8292$ and no information about location
or scale is provided; this is clearly inadequate.

We approach the problem via the method of *profile likelihood*
(see for example pg. 347 of @hinkley).
Initially, we treat the location and scale parameters as *nuisance parameters*
in order to obtain an estimate $\widehat{\nu}$ of the degrees-of-freedom parameter $\nu$.
With $\widehat{\nu}$ in hand, we then estimate location and scale using ordinary maximum likelihood.


# The basic model #
The data are assumed to be sampled from $Y = c + s\, X$, where $c$ and $s$ are location and scale parameters
respectively, and
where $X$ is t-distributed with $\nu$ degrees of freedom;
one writes $X \sim t_{\nu}$. Thus the pdf of $Y$ is given by:
$$ f_{Y}(y) = \frac{1}{s} f_{X} \left(\dfrac{y-c}{s}\right)\; \mbox{ where }\;
f_{X}(x) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sqrt{\nu \pi}\, \Gamma\left(\dfrac{\nu}{2}\right)}
\left(1 + \dfrac{x^2}{\nu}\right)^{-\dfrac{\nu+1}{2}}$$

For a given observation vector $y = [y_1, y_2,\ldots, y_m]^{\top}$ the likelihood function is then:

$$ l(\nu\,;\, c, s) = \sum_{i=1}^{m}\log \left(f_{X} \left(\dfrac{y_i-c}{s}\right)\right) - m\log(s) $$

The profile log-likelihood of $\nu$ is then

$$ l_{\mbox{prof}} (\nu) = \max_{c,\,s} {l(\nu\,;\,c,s)}  $$
and $\widehat{\nu}$ is obtained as:

$$ \widehat{\nu} = \arg\left(\max_{\nu}{l_{\mbox{prof}}(\nu)}\right) $$
and finally estimates $(\widehat{c},\, \widehat{s})$ are found as

$$ (\widehat{c},\, \widehat{s}) = \arg\left(\max_{c,\,s} l(\widehat{\nu}\,; c,s) \right) $$
The partial derivatives of $l(\nu\,; c, s)$ with respect to $c$ and $s$ will be used.

$$ \dfrac{\partial\, l(\nu\,;\, c, s)}{\partial\,c} =
\sum_{i=1}^{m} \dfrac{(\nu+1)\,(y_i-c)}
{(y_i -c)^2 + \nu\,s^2}
$$

$$ \dfrac{\partial\, l(\nu\,;\, c, s)}{\partial\,s}
= \sum_{i=1}^{m} \dfrac{\nu \left((c - y_i)^2 - s^2\right)}
{s\left((c-y_i)^2 + \nu\,s^2\right)}
$$

# Implementation
The above formulae are coded in R as follows:

```r
logpdf <- function (t, Location, Scale, nu) {
  temp = (t - Location) / Scale;
  temp = log (dt (temp, nu)) - log (Scale)
  temp
}

derv.Location.Scale <- function (t, Location, Scale, nu) {
  dwrts = (nu/Scale) * ((Location - t)^2 - Scale^2) / ((Location-t)^2 + nu * Scale^2)
  dwrtc = (nu + 1) * (t - Location) / ((t - Location)^2 + nu * Scale^2)
  c (dwrtc, dwrts)
}

loglk <- function (x, Location, Scale, nu) {
  temp = sum (logpdf (x, Location, Scale, nu))
  temp
}

score.Location.Scale <- function (x, Location, Scale, nu) {
  temp = derv.Location.Scale (x, Location, Scale, nu)
  temp = matrix (temp, nrow=length(x))
  res <- apply (temp, 2, sum)
  res 
}

target <- function (cands, nu, x) {
  0.0 - loglk (x=x, cands[1], cands[2], nu=nu)
}

gr <- function (cands, nu, x) {
  -1 * score.Location.Scale (x=x, cands[1], cands[2], nu=nu)
}

get.cands.fixednu <- function (nu, x, init_cands) {
  options (warn=-1)
  res = optim (par=init_cands, fn=target, gr=gr,
               method="BFGS",
               nu=nu,
               x=x)
  options (warn=0)
  res$par
}

nu.target <- function (nu, init_cands, x) {
  cands = get.cands.fixednu (nu, x, init_cands)
  target (cands = cands, nu=nu, x=x)
}
```
The next function `get.tparams.only` returns $\widehat{\nu}$, as well as $\widehat{c}$ and $\widehat{s}$,
but no confidence values yet; the latter will be provided later on.
Note that the optimization procedure requires initial guesses for $\nu$, $c$ and $s$.
The arguments `lower` and `upper` give a range in which the user believes $\nu$ lies, 
while the sample mean is used as the initial guess for $c$. Our initial guess
at the scale parameter is motivated by the discussion on page 86 of @serfling; we use
$$ s_{\mbox{init}} = \left(\widehat{F}_{m}^{-1}(3/4) - \widehat{F}_{m}^{-1}(1/4)\right) / 2 / 0.6745 $$
where $\widehat{F}_m^{-1}()$ is the empirical quantile function of the observation set $y$. 



```r
get.tparams.only <- function (x, lower, upper) {
  c_init = mean (x)
  s_init = 0.5 * (quantile (x, 0.75) - quantile (x, 0.25)) / 0.6745
  init_cands = c (c_init, s_init)
  names (init_cands) = NULL
  
  
  ans = optim (c((lower + upper)/2), fn=nu.target, method="Brent",
               lower = lower, upper = upper,
               init_cands = init_cands, x=x)
  
  tmp = get.cands.fixednu (ans$par, x, init_cands)
  ans$Location = tmp[1]
  ans$Scale = tmp[2]
  ans$nu = ans$par
  
  ans
}
```

Let us try `get.tparams.only` on our simulated data set $y$:

```r
my.estimate <- get.tparams.only (y, 2.1, 20)
str <- sprintf ("nu_hat = %.2f, c_hat = %.2f, s_hat = %.2f",
	my.estimate$nu, my.estimate$Location, my.estimate$Scale)
print (str)
```

```
## [1] "nu_hat = 4.93, c_hat = 1.81, s_hat = 2.94"
```
The true values are $5.3$, $2$ and $3$ respectively.

## Obtaining confidence intervals
Although in principle, confidence intervals for the MLE are available through the information
matrix, we found the particular form of $f_X(\cdot)$ too difficult to obtain second derivatives
for in closed form. Consequently, we resort to bootstrapping. 


```r
get.tparams.withCI  <- function(y, lower, upper, CL=0.95, B=99) {
  params <- get.tparams.only (y, lower, upper)
  bs.nu <- rep (0.0, B)
  bs.Location <- rep (0.0, B)
  bs.Scale <- rep(0.0, B)

  for (i in 1:B) {
	ind <- sample.int (n=length(y), size=length(y), replace=TRUE)
	bs.dat = y[ind]
	bs.ans = get.tparams.only (y[ind], lower=lower, upper=upper)
	bs.nu [i] = bs.ans$nu
	bs.Location [i] = bs.ans$Location
	bs.Scale [i] = bs.ans$Scale
  }

  params$boot.nu = mean (bs.nu)
  
  p = CL + (1-CL)/2
  m = mean (bs.nu)
  s = sd (bs.nu)
  params$nu.CI = params$par + c(-1,1) * quantile ((bs.nu - m)/s, p) * s / sqrt(B)
  
  m = mean (bs.Location)
  s = sd (bs.Location)
  params$Location.CI = params$Location + c(-1,1) * quantile ((bs.Location - m)/s, p) * s / sqrt(B)
  
  m = mean (bs.Scale)
  s = sd (bs.Scale)
  params$Scale.CI = params$Scale + c(-1,1) * quantile ((bs.Scale - m)/s, p) * s / sqrt(B)
  
  params
}
```
Now we test this on our simulated data $y$


```r
ptm <- proc.time()
CL <- 0.95
result <- get.tparams.withCI (y, 2.1, 20, CL = CL) 
proc.time() - ptm
```

```
##    user  system elapsed 
##  12.132   0.000  12.146
```

```r
str <- sprintf ("nu = %.2f, %.0f%% CI = [%.2f, %.2f]",
	result$nu, (100*CL), result$nu.CI[1], result$nu.CI[2])
print (str)
```

```
## [1] "nu = 4.93, 95% CI = [4.57, 5.29]"
```

```r
str <- sprintf ("Location = %.2f, %.0f%% CI = [%.2f, %.2f]",
	result$Location, (100*CL), result$Location.CI[1], result$Location.CI[2])
print (str)
```

```
## [1] "Location = 1.81, 95% CI = [1.77, 1.85]"
```

```r
str <- sprintf ("Scale = %.2f, %.0f%% CI = [%.2f, %.2f]",
	result$Scale, (100*CL), result$Scale.CI[1], result$Scale.CI[2])
print (str)
```

```
## [1] "Scale = 2.94, 95% CI = [2.91, 2.97]"
```


# References

