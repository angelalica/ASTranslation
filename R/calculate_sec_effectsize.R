#for a given SNP, returns the 7 values associated with secondary effect size
calculate_sec_effsize = function(RNA.mat, RNA.pat, RIBO.mat, RIBO.pat) {

  effsize.vec = vector(mode = "numeric", length = 7)

  sum.RNA.mat = sum(RNA.mat)
  sum.RNA.total = sum(RNA.mat, RNA.pat)

  sum.RIBO.mat = sum(RIBO.mat)
  sum.RIBO.total = sum(RIBO.mat, RIBO.pat)

  RNA.ratio = round((sum.RNA.mat / sum.RNA.total), 3)
  RIBO.ratio = round((sum.RIBO.mat / sum.RIBO.total), 3)

  diff = abs(RNA.ratio - RIBO.ratio)
  effsize.vec = c(diff, RNA.ratio, RIBO.ratio, sum.RNA.mat, sum.RNA.total, sum.RIBO.mat, sum.RIBO.total)
  return(effsize.vec)
}

# return a matrix where each row is the 7 secondary effect size values for each SNP
calculate_sec_effsize_all = function(RNA.all, RIBO.all, num.SNPs) {
  sec.effsize.all = matrix(nrow = num.SNPs, ncol = 7)
  colnames(sec.effsize.all) = c("Diff", "RNA.Ratio", "RIBO.Ratio", "mat.RNA", "all.RNA", "mat.RIBO", "all.RIBO")

  for (i in 1:num.SNPs) {
    sec.effsize.all[i, ] = calculate_sec_effsize(RNA.all[2*i-1, ], RNA.all[2*i, ], RIBO.all[2*i-1, ], RIBO.all[2*i, ])
  }

  return(sec.effsize.all)
}
