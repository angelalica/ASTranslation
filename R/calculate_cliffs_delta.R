#for a given SNP, calculate Cliff's Delta between the RNA and RIBO maternal / total ratio distributions
calculate_cliffs_delta = function(RNA.ratio, RIBO.ratio, num.sims) {
  combined = c(RNA.ratio, RIBO.ratio)
  c.delta = cliff.delta(combined, c(rep("RNA", num.sims), rep("RIBO", num.sims)))
  result = c.delta$estimate
  return(round(abs(result), 3))
}

# return a vector of Cliff Delta values for each SNP
calculate_cliffs_delta_all = function(RNA.counts.ratios, RIBO.counts.ratios, num.SNPs, num.sims) {

  c.delta.all = vector(mode = "numeric", length = num.SNPs)
  for (i in 1:num.SNPs) {
    c.delta.all[i] = calculate_cliffs_delta(RNA.counts.ratios[i, ], RIBO.counts.ratios[i, ], num.sims)
  }
  return(c.delta.all)
}
