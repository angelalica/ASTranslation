#returns all bootstrapped Ch1+/Ch2+ ratios for one fraction
bootstrap_ratio = function (x, nsims) {
  observation = c(rep(1, x[1]), rep(2, x[2]), rep(0, x[3]) )
  ratios = c()
  for (j in 1:nsims) {
    boot_tmp = sample(observation, replace =T)
    ratios = c(ratios , length(which(boot_tmp == 1)) / length(which(boot_tmp == 2))  )
  }
  return (ratios)
}
