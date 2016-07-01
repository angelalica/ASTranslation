# bootstrap ratios each fraction, generates p-values, returns all bootstrapped Ch1+/Ch2+ ratios for each fraction
significance_testing = function (condensed.data, expected.value, n, zoom.range, output.name) {

  ratio.matrix = NULL
  all.pvalues = NULL

  #generate pvalues for each fraction from the bootstrap values
  for (i in 1:dim(condensed.data)[1]) {
    ratios = bootstrap_ratio(condensed.data[i, 1:3], 1000)
    ratio.matrix = cbind(ratio.matrix, ratios)

    num.extreme = ratios[which(ratios <= expected.value)]
    pvalue = length(num.extreme) / 1000 #calculate the p-values
    if (pvalue > 0.5) { pvalue = 1 - pvalue }
    pvalue = pvalue * 2 #double sided test

    all.pvalues = append(all.pvalues, pvalue)
  }

  #prepare final output
  final.output = cbind(all.pvalues, round(condensed.data[, 5], digits = 3))
  colnames (final.output) = c("p-values", "Proportion of Total Concentration")
  write.csv(final.output, output.name, row.names = rownames(condensed.data))

  colnames(ratio.matrix) = rownames(condensed.data)

  # if there is a specific zoom, return only the relevant fractions for plotting
  if (!is.null(zoom.range)) { return (ratio.matrix[, zoom.range[1]:zoom.range[2]])}
  return (ratio.matrix)
}
