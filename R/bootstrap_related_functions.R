#creates a bootstrap distribution of maternal / total RNA or RIBO counts ratio for one SNP by
#sampling the counts for each allele and summing them up to form the maternal / total ratio
#input: counts vector (either RNA or RIBO) from maternal and paternal allele,
#psuedocounts, and number of times to bootstrap
#output: vector of maternal / total counts ratio for one SNP (either RNA or RIBO)
bootstrap = function (mat.vec, pat.vec, pcounts, nsims) {
  row = vector(length = nsims)
  for (i in 1:nsims) {
    #sample the appropriate positions first
    sample.positions = sample.int(length(mat.vec), size = length(mat.vec), replace = TRUE)

    #apply the same positions to both mat and pat allele
    mat.samplevec = mat.vec[c(sample.positions)]
    pat.samplevec = pat.vec[c(sample.positions)]
    row[i] = (sum(mat.samplevec) + pcounts) / (sum(mat.samplevec) + sum(pat.samplevec) + 2*pcounts)
  }
  return (row)
}

#bootstrap all SNPs
#returns a matrix where each row is all of the bootstrapped ratios generated for a particular SNP
generate_bootstrapped_ratios = function(RNA.all, RIBO.all, num.SNPs, nsims, pcounts) {

  RNA.counts.ratios = matrix(nrow = num.SNPs, ncol = nsims)
  RIBO.counts.ratios = matrix(nrow = num.SNPs, ncol = nsims)

  #bootstrap each SNP
  for (i in 1:num.SNPs) {
    RNA.counts.ratios[i, ] = bootstrap(RNA.all[(2*i-1), ], RNA.all[(2*i), ], pcounts, nsims)
    RIBO.counts.ratios[i, ] = bootstrap(RIBO.all[(2*i-1), ], RIBO.all[(2*i), ], pcounts, nsims)
  }
  return(list(RNA.counts.ratios, RIBO.counts.ratios))
}

#print out warnings about if median of bootstrapped ratios for any SNP is suspiciously close to 0 or 1
warn_extreme_ratios = function(RNA.counts.ratios, RIBO.counts.ratios, num.SNPs, names.SNPs) {
  for (i in 1:num.SNPs) {
    m.RNA = median(RNA.counts.ratios[i, ])
    m.RIBO = median(RIBO.counts.ratios[i, ])
    if (m.RNA < 0.1 || m.RIBO < 0.1) {
      cat(sprintf("Warning: mean bootstrapped maternal/total ratio counts for either RNA or RIBO in %s is less than 0.1\n", names.SNPs[i]))
    }
    if (m.RNA > 0.9 || m.RIBO > 0.9) {
      cat(sprintf("Warning: mean bootstrapped maternal/total ratio counts for either RNA or RIBO in %s is greater than 0.9\n", names.SNPs[i]))
    }
  }
}
