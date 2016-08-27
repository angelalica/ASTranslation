#removes the SNPs with zero reads from both the RNA and RIBO counts files
#returns new (filtered) RNA.all and RIBO.all (as a list of two matrices)

filter_zero_read_SNPs = function(RNA.all, RIBO.all) {

  SNPs.RNA.zero = identify_zero_read_SNPs(RNA.all, "RNA-seq", 12) #12 = length of string "*aternalRNA"
  SNPs.RIBO.zero = identify_zero_read_SNPs(RIBO.all, "ribosome profiling", 13)

  #get all SNPs that have missing counts for either data table
  SNPs.to.remove = union(SNPs.RNA.zero, SNPs.RIBO.zero)

  #actually remove these SNPs from consideration
  SNPs.all = substr(rownames(RNA.all), 1, nchar(rownames(RNA.all)) - 12)
  indices.to.remove = which(!is.na(match(SNPs.all, SNPs.to.remove))) #match
  RNA.all = RNA.all[-(indices.to.remove),]
  RIBO.all = RIBO.all[-(indices.to.remove),]

  return(list(RNA.all, RIBO.all))
}

#identify the names of the SNPs that have zero reads
identify_zero_read_SNPs = function(raw.counts.matrix, tech, num.remove) {
  counts.sum = rowSums(raw.counts.matrix)
  indices = which(counts.sum == 0)

  cat(sprintf("\nThe following SNP alleles have zero %s reads:\n", tech))
  cat(names(indices), sep = "\n")

  #get rid of the suffix (i.e. ".maternalRNA" so that just the SNP name remains)
  SNPs.zero = substr(names(indices), 1, nchar(names(indices)) - num.remove)
  return(unique(SNPs.zero))
}

#creates a bootstrap distribution of maternal / total RNA or RIBO counts ratio for one SNP by
#sampling the counts for each allele and summing them up to form the maternal / total ratio
#input: counts vector (either RNA or RIBO) from maternal and paternal allele,
#       psuedocounts, and number of times to bootstrap
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
