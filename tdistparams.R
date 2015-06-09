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
