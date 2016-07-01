#get overall Ch1+/Ch2+ ratio (expected value)
calculate_totals = function(condensed.data) {
  totals.calc = apply(condensed.data, 2, function(x) {sum(x)}) #sum up both, mutant, and wildtype columns
  total.ratio.calc = totals.calc[1] / totals.calc[2]
  return(total.ratio.calc)
}

#group data by ribosomal fractions (number of ribosomes)
group_data = function(condensed.data, grouped.fractions) {

  grouped.data = NULL

  #average the experimental fractions that correspond to each ribosomal fraction
  for (i in 1:length(grouped.fractions)) {
    if (grouped.fractions[i] > 1) {
      temp = apply(as.matrix(condensed.data)[1:grouped.fractions[i], 1:4], 2, FUN = mean)
      temp.conc = sum(condensed.data[1:grouped.fractions[i], 5])
      temp = append(temp, temp.conc)
    } else { temp = condensed.data[1, ] }
    condensed.data = condensed.data[-(1:grouped.fractions[i]), ]
    grouped.data = rbind(grouped.data, temp)
  }
  rownames(grouped.data) = c(0:(length(grouped.fractions) - 1))

  return(grouped.data)
}
