#gives the parameters for sizing and binning the overlapping histogram appropriately
find_histogram_parameters = function(RNA.ratio, RIBO.ratio) {
  info.RNA = hist(RNA.ratio, plot = FALSE)
  info.RIBO = hist(RIBO.ratio, plot = FALSE)

  min.x = min(info.RNA$breaks[1], info.RIBO$breaks[1])
  max.x = max(info.RNA$breaks[-1], info.RIBO$breaks[-1])
  width = min(info.RNA$breaks[2] - info.RNA$breaks[1], info.RIBO$breaks[2] - info.RIBO$breaks[1])
  len = ceiling((max.x-min.x) / width)

  new.hist.RNA = hist(RNA.ratio, breaks = seq(min.x, max.x, length = len), plot = FALSE)
  new.hist.RIBO = hist(RIBO.ratio, breaks = seq(min.x, max.x, length = len), plot = FALSE)
  y.max = max(max(new.hist.RNA$counts), max(new.hist.RIBO$counts))

  result = c(min.x, max.x, len, y.max)
  return (result)
}

#makes the histogram plot for one particular SNP
plot_one_histogram = function(RNA.ratio, RIBO.ratio, c.delta, p.value, diff, SNP) {
  params = find_histogram_parameters(RNA.ratio, RIBO.ratio)
  sub = sprintf("C.Delta = %s, p.value = %s, Diff = %s", c.delta, p.value, diff)

  hist(RNA.ratio, breaks = seq(params[1], params[2], length = params[3]),
       ylim = c(0, params[4]), col=rgb(1,0,0,0.5),
       main= sprintf("%s \n %s", SNP, sub), cex.main = 1.1,
       xlab="Maternal/Total Counts")

  hist(RIBO.ratio, breaks = seq(params[1], params[2], length = params[3]),
       col=rgb(0,0,1,0.5), add = T)

  legend("topright", col=c(rgb(1,0,0,0.5), rgb(0, 0, 1, 0.5)),
         lty = c(1, 1), legend = c("RNA","Ribo"), cex = 0.5,
         inset = c(-0.05, 0), xpd = TRUE)
}

#plots histograms of each candidate SNP and outputs all of the plots to one pdf file
plot_all_histograms = function(RNA.counts.ratios, RIBO.counts.ratios, final.candidates, ordered.indices, output.pdf) {
  pdf (file = output.pdf)
  par(mfrow = c(2, 2))
  j = 1
  for (i in ordered.indices) {
    p.value = round(final.candidates[j, 1], 3)
    c.delta = final.candidates[j, 2]
    diff = final.candidates[j, 3]
    SNP = rownames(final.candidates)[j]

    plot_one_histogram(RNA.counts.ratios[i, ], RIBO.counts.ratios[i, ],
                       c.delta, p.value, diff, SNP)
    j = j+1
  }
  dev.off()
}
