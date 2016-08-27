#' Determines whether a candidate SNP has allele specific translation.
#'
#' This function helps the user determine whether a particular heterozygous SNP
#' has allele-specific translation based on data from an assay that combines
#' polysome profiling and digital droplet PCR (ddPCR). It uses bootstrapping
#' techniques to determine whether the amount of mutant and wildtype transcripts
#' in each fraction (directly from the polysome profiling or grouped by number of
#' ribosomes) is significantly different from the expected values. Each fraction
#' corresponds to a particular weight (heavier transcripts contain a greater number
#' of ribosomes, indicating greater translation).
#'
#' @param file.name Filename of the ddPCR counts data file to be read in. The
#'   data must be formatted in the following manner: each row corresponds to a
#'   ddPCR well, and there must be exactly 7 columns - well ID, fraction ID,
#'   Ch1+/Ch2+ counts, Ch1+/Ch2- counts, Ch1-/Ch2+ counts, Ch1-/Ch2- counts, and
#'   proportion of total concentration.
#' @param gene.name A string indicating the name of gene being studied.
#' @param expected.value A number indicating the overall Ch1+/Ch2+ ratio that
#'   would be expected if there were no allele specific translation. The user
#'   can include if this value is already known. If \code{NULL}, the function will
#'   automatically calculate an expected value by summing over all fractions.
#'   Defaults to \code{NULL}.
#' @param grouped A boolean indicating whether to group fractions by number of
#'   polysomes. Defaults to \code{FALSE}.
#' @param grouped.fractions A vector indicating the number of experimental
#'   fractions corresponding to each ribosomal fraction. Defaults to \code{NULL}.
#' @param zoom.range A vector of 2 numbers that indicate the particular range of
#'   experimental fractions to graph. Both numbers must be positive. The first number
#'   must be strictly smaller than the second number and neither can be greater than the
#'   total number of fractions. Defaults to \code{NULL}, which corresponds to
#'   plotting all fractions. NOTE: zoom.range can only be used with fractions that are
#'   NOT grouped by ribosomes.
#' @param nsims The number of bootstrap simulations to run for each fraction
#' @param output.pvals File name of the csv output file
#' @param output.plots File name of the pdf output file
#' @return There are two outputs. The first is a CSV file with each fraction /
#'   group, it's p-value, and its proportion of the total concentration. The
#'   second output is a PDF containing two aligned graphs. One graph corresponds
#'   to a boxplot of all bootstrapped Ch1+/Ch2+ ratios relative to the expected
#'   value (either calculated or inputted by the users). This makes it visually
#'   clear which fractions are significantly different (p < 0.05 for a
#'   two-tailed test) from the expected ratios. The other graph is a barplot
#'   that graphs the ddPCR concentrations from each fraction so that users can
#'   see the distribution and relative concentrations of RNA across all
#'   fractions (from the early fractions with no ribosomes to the heavier
#'   fractions with 5+ ribosomes).
#' @examples
#'   test_candidates (filename = "gene1_ddPCR_data.csv", gene.name = "GENE1")
#'   test_candidates (filename = "gene2_ddPCR_data.csv", gene.name = "GENE2", zoom.range = c(5, 16))
#'   test_candidates (filename = "gene3_ddPCR_data.csv", gene.name = "GENE3", expected.value = 1.2, grouped = TRUE, grouped.fractions = c(7, 3, 1, 1, 1, 5))

# overall function (that user interacts with), provides significance testing for each fraction and generates summary plots
test_candidates = function (file.name,
                    gene.name,
                    expected.value = NULL,
                    grouped = FALSE,
                    grouped.fractions = NULL,
                    zoom.range = NULL,
                    nsims = 1000,
                    output.pvals = "output_pvals.csv",
                    output.plots = "output_plots.pdf") {

  #checking all arguments first
  if(!is.character(file.name)) {stop ("file.name must be a string")}
  if(!is.character(gene.name)) {stop ("gene.name must be a string")}

  if(grouped && is.null(grouped.fractions)) {stop ("Please enter a vector of numbers for grouped.fractions")}
  if(!grouped && !is.null(grouped.fractions)) {stop ("Please set grouped = TRUE")}

  if(grouped && !is.null(zoom.range)) {stop ("Cannot zoom when fractions are grouped by ribosomes")}
  if((zoom.range[1] <= 0) || (zoom.range[2] <= 0)) {stop ("zoom.range values cannot be less than or equal to 0")}
  if(zoom.range[1] >= zoom.range[2]) {stop ("Please enter an appropriate range (first value must be strictly smaller than second value).")}
  if(!is.character(output.pvals) || substr(output.pvals, nchar(output.pvals)-3, nchar(output.pvals)) != ".csv") {
    stop("output.pvals must be a string and a valid .csv file")
  }
  if(!is.character(output.plots) || substr(output.plots, nchar(output.pvals)-3, nchar(output.pvals)) != ".pdf") {
    stop("output.plots must be a string and a valid .csv file")
  }

  #condensed_data has 5 columns: Ch1+, Ch2+, Ch1-/Ch2-, Ch1+/Ch2+ ratio, Proportion of Total Concentration
  condensed.data = process_data(file.name)

  #which expected value to use
  if (!is.null(expected.value)) {
    #if user has a pre-calculated expected value
    expected = expected.value
  } else {
    #if user wants to sum up totals to determine an expected value
    expected = calculate_totals(condensed.data[, 1:3])
  }

  #if grouped by ribosomal fractions
  if (grouped) {
    condensed.data = group_data(condensed.data, grouped.fractions)
  }

  if (zoom.range[1] > nrow(condensed.data) || zoom.range[2] > nrow(condensed.data)) {
    stop("zoom.range values are out of bounds")
  }
  #do bootstrapping for each fraction and return the all bootstrapped ratios for each fraction
  ratio.matrix = significance_testing(condensed.data, expected, nsims, zoom.range, output.pvals)

  # plot the two output graphs
  if (!is.null(zoom.range)) {
    #if user wants a specific zoom
    plot_ratios_and_proportions(condensed.data[zoom.range[1]:zoom.range[2], 5], ratio.matrix, expected, grouped, gene.name, output.plots)
  } else {
    #plot everything
    plot_ratios_and_proportions(condensed.data[, 5], ratio.matrix, expected, grouped, gene.name, output.plots)
  }
}

