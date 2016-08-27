#screen each SNP based on whether Cliffs' Delta and secondary effect size pass user-defined threshold
#returns a matrix where each row is a SNP that passed the screen, along with its associated
#Cliff's delta and secondary effect size values
screen_candidates = function(cliffs.all, sec.effsize.all, names.SNPs, c.threshold, d.threshold) {
  effsize.all = cbind(cliffs.all, sec.effsize.all)
  rownames(effsize.all) = names.SNPs
  effsize.candidates = effsize.all[which(effsize.all[, 1] > c.threshold), ]
  effsize.candidates = effsize.candidates[which(effsize.candidates[, 2] > d.threshold), ]

  return(effsize.candidates[order(effsize.candidates[, 1], decreasing = TRUE), ])
}
