#' Determining which SNPs are candidates for allele specific translation
#' This function helps users identify which genetic variants are candidates for
#' allele-specific translation (i.e. may lead to differential translation) by
#' using RNA-Seq and Ribosome profiling data from all heterozygous SNPs in the
#' transcriptome of a given individual.
#'
#' @param RNA.file Filename of the RNA-Seq counts for each SNP. Each row is either the maternal
#'        or paternal allele of a SNP and its associated number of RNA-Seq reads at each position.
#'        The number of columns is defined by the size of the window (i.e. for RNA, usually around
#'        75 nucleotides.
#'
#' @param RIBO.file Filename of the Ribosome profiling counts for each SNP. See \code{RNA.file}
#'        for more details.
#'
#' @param output.csv Filename of the output csv file where the numerical output (effect sizes
#'        and p values) will be stored.
#'
#' @param output.pdf Filename of the output pdf file where the graphical output (histograms
#'        of the RNA-Seq and ribosome profiling bootstrapped ratio distributions) will be
#'        stored
#'
#' @param alpha Threshold for what p-values are considered significant. Default is 0.05.
#'
#' @param pcounts Pseudocounts for regularization when calculating the maternal / total ratio
#'        counts (for each bootstrapped sample) so tha the RNA-Seq ratio and ribosomal
#'        profiling ratio can be comparable. Default value is 50. Should be a value between
#'        ~30 (the window size ribosome profiling) and ~75 (the window size for RNA-seq).
#'
#' @param c.threshold Threshold for what Cliff's Delta values are considered significant.
#'        Default is 0.9.
#'
#' @param d.threshold Threshold for what secondary effect size difference is significant.
#'        Default is 0.1 (i.e. a minimum of a 10% difference between the allelic ratio
#'        of the RNA-seq counts and ribosome profiling counts).
#'
#' @param num.sims Number of times to bootstrap. Default is 1000.
#'
#' @return Two files. 1) A .csv file with all candidate SNPs and their associated
#'        p-values, Cliff's Delta, and secondary effect size values (9 columns total).
#'        2) A .pdf file with a histogram of the bootstrapped ratio distributions
#'        for RNA-Seq (red) and ribosome profiling (blue) for each candidate SNP. This makes
#'        it possible to easily visualize how much overlap there are between the
#'        two distributions.
#'
#' @examples find_candidates("RNA.all.csv", "RIBO.all.csv", "final.effsize.csv", "final.plots.pdf")
#' find_candidates("RNA.all.csv", "RIBO.all.csv", "final.effsize.csv", "final.plots.pdf", num.sims = 5000, c.threshold = 0.5)
#'
find_candidates = function(RNA.file,
                           RIBO.file,
                           output.csv,
                           output.pdf,
                           alpha = 0.05,
                           pcounts = 50,
                           c.threshold = 0.9,
                           d.threshold = 0.1,
                           num.sims = 1000) {

  filtered.data = read_and_filter(RNA.file, RIBO.file)
  RNA.all = filtered.data[[1]]
  RIBO.all = filtered.data[[2]]

  num.SNPs = dim(RNA.all)[1] / 2
  names.SNPs = unique(substr(rownames(RNA.all), 1, nchar(rownames(RNA.all)) - 12))

  #bootstrap the counts for RNA and RIBO to create a distribution of allelic ratios (maternal / total)
  ratios.counts = generate_bootstrapped_ratios(RNA.all, RIBO.all, num.SNPs, num.sims, pcounts)
  RNA.counts.ratios = ratios.counts[[1]]
  RIBO.counts.ratios = ratios.counts[[2]]

  warn_extreme_ratios(RNA.counts.ratios, RIBO.counts.ratios, num.SNPs, names.SNPs)

  #perform 2 different calculations of effect size
  cliffs.all = calculate_cliffs_delta_all(RNA.counts.ratios, RIBO.counts.ratios, num.SNPs, num.sims)
  sec.effsize.all = calculate_sec_effsize_all(RNA.all, RIBO.all, num.SNPs)
  possible.candidates = screen_candidates(cliffs.all, sec.effsize.all, names.SNPs, c.threshold, d.threshold)

  #get pvalues of screened candidates
  ordered.indices = match(rownames(possible.candidates), names.SNPs)
  p.values = calculate_pvalue_all(RNA.counts.ratios, RIBO.counts.ratios, ordered.indices)
  final.candidates = screen_pvalues(p.values, possible.candidates, alpha)

  #return results (table of significant SNPs with p-values and all effect sizes reported for each)
  write.csv(final.candidates, output.csv)

  final.ordered.indices = match(rownames(final.candidates), names.SNPs)
  plot_all_histograms(RNA.counts.ratios, RIBO.counts.ratios, final.candidates, final.ordered.indices, output.pdf)
}
