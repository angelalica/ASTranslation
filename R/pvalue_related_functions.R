#for a given SNP, calculates p-value based on the largest confidence level that would give rise to
#non-overlapping confidence intervals between the RNA and RIBO bootstrapped ratio distributions
calculate_pvalue = function (RNA.ratio, RIBO.ratio) {

  i = 99 #start at 99 percent confidence interval
  while(TRUE) {
    alpha = 100 - i
    l.percentile = (alpha/2) / 100
    h.percentile = (100 - alpha/2) / 100

    RNA.lower = quantile(RNA.ratio, l.percentile)
    RNA.upper = quantile(RNA.ratio, h.percentile)

    RIBO.lower = quantile(RIBO.ratio, l.percentile)
    RIBO.upper = quantile(RIBO.ratio, h.percentile)

    if (mean(RIBO.ratio) > mean(RNA.ratio)) {
      if (RIBO.lower > RNA.upper) { break }
    } else {
      if (RNA.lower > RIBO.upper) { break }
    }
    i = i - 1
  }
  return((100 - i) / 100)
}

#return a vector of p-values for all SNPs that pass the thresholds (candidates)
calculate_pvalue_all = function (RNA.counts.ratios, RIBO.counts.ratios, ordered.indices) {
  p.values = vector(mode = "numeric", length = length(ordered.indices))
  j = 1
  for (i in ordered.indices) {
    p.values[j] = calculate_pvalue(RNA.counts.ratios[i, ], RIBO.counts.ratios[i, ])
    j = j + 1
  }
  return(p.values)
}

#return a table of final candidates with p-values and effect sizes (sorted based on pvalue)
screen_pvalues = function(p.values, possible.candidates, alpha) {
  possible.candidates = cbind(p.values, possible.candidates)
  final.candidates = possible.candidates[which(possible.candidates[, 1] <= alpha), ]
  return(final.candidates[order(final.candidates[, 1]), ])
}
