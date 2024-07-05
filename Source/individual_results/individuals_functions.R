

## CREATE PER SAMPLE FILTERED LIST OF ANNOTATED REPEATS
by_sample = function(reads, savename, reference, min_sd=3, min_length = 50){
  # read in (filtered) experimental data:
  df = read.csv(file=reads, sep="\t", header = TRUE)
  # read in precompute reference data
  ref = read.csv(file=reference, sep="\t", header = TRUE)
  # merge both
  df <- merge(df, ref[c("Call_ID", "Med_Units", "Mean_Units", "SD","Instability_Rating","Status", "Gene", "Region")], by="Call_ID")
  write.table(df, paste0(output, "_all_samples_full.csv"), sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # TODO : USE A PRECALCULATED REFSET TO GET SD/MEAN VALUES
  # select rows i df WHERE either : 
  #   - df$All1 > min(min_length, df$Mean_Units + min_sd * df$SD)
  #   - df$All2 > min(min_length, df$Mean_Units + min_sd * df$SD)
  # assign a large SD to SD == NA
  df$SD[is.na(df$SD)] <- 100000
  df_sign = subset(df, df$All1 > min(min_length, df$Mean_Units + min_sd * df$SD) | df$All2 > min(min_length, df$Mean_Units + min_sd * df$SD))
  write.table(df_sign, paste0(output, "_all_samples_hard_cutoff.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
  # create output files per sample: 
  for (sample in unique(df$Sample_ID)){
    df_sample = subset(df_sign, df_sign$Sample_ID == sample)
    write.table(df_sample, paste0(output, "_", sample, ".csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
  }

}
