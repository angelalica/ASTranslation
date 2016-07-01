#cleans up and processes input data into correct form for downstream analysis
#issues warnings to the user if the data has possible calling errors
process_data = function(file.name) {

  #read in original data table
  full.data = read.table(file.name, header = T, sep = ",")

  #warn user if double positives in a particular well are greater than 90% of total droplets
  counts.data = unique(full.data[, 1:6])
  num.fractions = length(rownames(counts.data))
  for (i in 1:num.fractions) {
    if ((counts.data$Ch1.Ch2.[i] / sum(counts.data[i, 3:6])) > 0.9) {
      cat(sprintf("WARNING: Double positive droplets in Well %s (%s) are >90 percent of total droplets\n", counts.data$Well[i], counts.data$X[i]))
    }
  }

  #add Ch1+/Ch2+ droplets to both Ch1+ and Ch2+
  for (i in 1:length(rownames(full.data))) {
    full.data[i,4] = full.data[i,4] + full.data[i,3]
    full.data[i,5] = full.data[i,5] + full.data[i,3]
    full.data[i,8] = c(full.data[i,4] / full.data[i,5])
  }

  #get rid of the columns corresponding to each well (no longer needed)
  full.data = full.data[, -1]
  colnames(full.data) = c("fraction", "Ch1+/Ch2+", "Ch1+/Ch2-",
                          "Ch2+/Ch1-", "Ch1-/Ch2-",
                          "Proportion of Total Concentration", "MT/WT ratio")

  #initialize empty data frame for the output / processed data
  condensed.data = NULL

  #handle the replicates through averaging and summing
  for (fraction in unique(full.data[, 1])) {
    condensed.row = c(mean(full.data[which(full.data[, 1] == fraction), 3]), #Ch1+/Ch2-
                      mean(full.data[which(full.data[, 1] == fraction), 4]), #Ch2+/Ch1-
                      mean(full.data[which(full.data[, 1] == fraction), 5]), #Ch1-/Ch2-
                      mean(full.data[which(full.data[, 1] == fraction), 7]), #Ch1+/Ch2+ Ratio
                      sum(full.data[which(full.data[, 1] == fraction), 6]))  #Conc Proportion
    rbind(condensed.data, condensed.row) -> condensed.data
  }

  rownames(condensed.data) = as.character(unique(full.data[, 1]))
  colnames(condensed.data) = c("Ch1", "Ch2", "neither", "Ch1+/Ch2+ ratio", "Fractional Concentration")

  return(condensed.data)
}
