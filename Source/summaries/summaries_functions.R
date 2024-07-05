by_locus = function(df, output, savename, build, humandb, globdir){
  #     Housekeeping
  df$All1 <- ifelse(df$All1 == 0, df$All2, df$All1)
  df$All2 <- ifelse(df$All2 == 0, df$All1, df$All2)
  df$EE50 <- ifelse(df$All1 >= 50 | df$All2 >= 50, 1, 0)
  df$EE100 <- ifelse(df$All1 >= 100 | df$All2 >= 100, 1, 0)
  df$avealleles <- (df$All1+df$All2)/2
  
  #     Summarize
  mdata <- melt(df[c("Call_ID", "All1", "All2")], id="Call_ID")
  medi <- aggregate(list("Median" = mdata$value), by=list("Call_ID" = mdata$Call_ID), median)
  medi <- medi[order(medi$Call_ID),]
  df <- merge(df,medi,by="Call_ID")
  df$stab <- ifelse(df$All1!=df$Median & df$All2!=df$Median, 1,
                    ifelse(df$All1==df$Median & df$All2!=df$Median, 0.5,
                           ifelse(df$All1!=df$Median & df$All2==df$Median, 0.5, 0)))
  ref <- aggregate(list("Ref_Units"=df$Ref_Units), by = list("Call_ID"=df$Call_ID, "Chr"=df$Chr, "Start"=df$Start), FUN = mean)
  avera <- aggregate(list("Mean"=df$avealleles), by = list("Call_ID"=df$Call_ID), FUN = mean)
  avera$Mean <- round(avera$Mean)
  counts <- aggregate(list("Count"=df$Call_ID), by = list("Call_ID"=df$Call_ID), FUN = length)
  stability <- aggregate(list("Stab"=df$stab), by = list("Call_ID"=df$Call_ID), FUN = mean)
  maxreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("Call_ID"=df$Call_ID), FUN = max)
  maxreps$max <- pmax(maxreps$All1, maxreps$All2)
  minreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("Call_ID"=df$Call_ID), FUN = min)
  minreps$min <- pmin(minreps$All1, minreps$All2)
  stddev <-  aggregate((All1+All2)~Call_ID, df, sd)
  ee50 <- aggregate(list("EE50"=df$EE50), by = list("Call_ID"=df$Call_ID), FUN = sum)
  ee100 <- aggregate(list("EE100"=df$EE100), by = list("Call_ID"=df$Call_ID), FUN = sum)
  rm(df)
  compile <- Reduce(function(x, y) merge(x, y, all=TRUE), list(ref, avera, medi, minreps[c("Call_ID", "min")], maxreps[c("Call_ID", "max")], counts, stability, stddev, ee50, ee100))
  compile$Status <- ifelse(compile$Stab == 0, "STABLE", "POLYMORPHIC")
  #     Annotate Repeats
  anno_dir <- paste0(output, "/Annovar")
  if (!dir.exists(anno_dir)) {
    system(paste0("mkdir ", anno_dir))
  }
  anno_out <- paste0(anno_dir, "/", savename)
  anno <- compile[c("Chr", "Start")]
  anno$End <- (anno$Start)+(3*compile$Ref_Units-1)
  anno$RA <- 0
  anno$AA <- "-"
  write.table(anno, paste0(anno_out, "_Anno.csv"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  # mk command explicit.
  perl_script = file.path(globdir,"Source", "summaries", "annotate_variation.pl")
  outfile_prefix = file.path(anno_dir,paste0(savename,"_Annovar"))
  outfile_name = paste0(anno_out,"_Anno.csv")
  cmd = paste("perl", perl_script, "--geneanno", "--outfile", outfile_prefix, "--buildver", paste0("hg", build), outfile_name, humandb)
  system(cmd)
  genes <- read.csv(paste0(anno_out, "_Annovar.variant_function"), sep = '\t', header = FALSE)
  genes$V1 <- revalue(genes$V1, c("upstream;downstream"="intergenic", "splicing"="intronic", "ncRNA_exonic"="ncRNA", "ncRNA_intronic"="ncRNA"))
  genes$V2 <- gsub("\\s*\\([^\\)]+\\)", "", genes$V2)
  genes$V2 <- gsub(",", "->", genes$V2)
  genes$V2 <- gsub(";", "->", genes$V2)
  genes$V7 <- paste(genes$V3, genes$V4, sep=".")
  genes$V7 <- gsub("chr", "", genes$V7)
  genes <- genes[c("V7", "V2", "V1")]
  names(genes) <- c("Call_ID", "Gene", "Region")
  compile <- merge(compile, genes, by="Call_ID")
  names(compile) <- c("Call_ID", "Chr", "Start", "Ref_Units", "Mean_Units", "Med_Units", "Min_Units", "Max_Units", "Hits", "Instability_Rating", "SD", "EE50", "EE100", "Status", "Gene", "Region")
  compile$Instability_Rating <- round(compile$Instability_Rating, 3)
  compile$SD <- round(compile$SD, 3)
  write.table(compile, paste0(output, savename, "_by_locus.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
  return(compile)
}


ratio_of_stability = function(df, output) {
  df <- df[c("Call_ID", "Mean_Units", "Med_Units", "Instability_Rating", "Status")]
  df$Med_Units <- round(df$Med_Units)
  unstable_reps <- subset(df, df$Status != "STABLE")
  stable_reps <- subset(df, df$Status == "STABLE")
  if (nrow(stable_reps) == 0) {
    count_stable <- data.frame(Med_Units = numeric(), Stable_Reps = numeric())
  } else {
    count_stable <- aggregate(list("Stable_Reps"=stable_reps$Med_Units), by = list("Med_Units"=stable_reps$Med_Units), FUN = length)
  } 
  if (nrow(unstable_reps) == 0) {
    count_unstable <- data.frame(Med_Units = numeric(), Unstable_Reps = numeric())
  } else {
    count_unstable <- aggregate(list("Unstable_Reps"=unstable_reps$Mean_Units), by = list("Med_Units"=unstable_reps$Med_Units), FUN = length)
  }
  
  count_total <- merge(count_stable, count_unstable, all = TRUE)
  count_total[is.na(count_total)] <- 0
  count_total$Total_Reps <- count_total$Stable_Reps+count_total$Unstable_Reps
  count_total$Flagged_Unstable <- (count_total$Unstable_Reps)/count_total$Total_Reps
  
  write.table(count_total, paste0(output, "_by_heterozygosity.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}

by_sample = function(df, locus, output){
  # POPULATION BASED
  df <- merge(df, locus[c("Call_ID", "Med_Units", "Mean_Units", "SD","Instability_Rating","Status", "Gene", "Region")], by="Call_ID")  
  poly <- subset(df, df$All1 != Med_Units | df$All2 != Med_Units)
  pa <- aggregate(poly$Sample_ID, by = list(poly$Sample_ID), FUN = length)
  a <- aggregate(df$Sample_ID, by = list(df$Sample_ID), FUN = length)
  b <- aggregate((All1+All2)/2 ~Sample_ID, data = df, FUN = mean)
  mdata <- melt(df[c("Sample_ID", "All1", "All2")], id="Sample_ID")
  c <- aggregate(mdata$value, by = list(mdata$Sample_ID), FUN = max)
  z <- data.frame("Sample_ID"=a$Group.1, "Highest_Repeat"=c$x, "Total_Repeats"=a$x, "Mean_Rep_length"=round(b$"(All1 + All2)/2", 2), "Unstable_Reps"=pa$x, "Percentage_Unstable_Reps"=round((pa$x/a$x*100), 2))
  z <- z[order(z$Sample_ID),]
  write.table(z, paste0(output, "_by_sample.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

}

by_chr = function(df, locus, output){
  # TODO : use a genome dict for this
  if (grepl("chr", as.character(df$Chr[1]), fixed = TRUE) == F){chrmbp <- data.frame("Chr"=c( as.character(c(1:22)),"X","Y"), "Mbp"=c(249,237,192,183,174,165,153,135,132,132,132,123,108,105,99,84,81,75,69,63,54,57,141,60))}
  else if (grepl("chr", as.character(df$Chr[1]), fixed = TRUE) == T) {chrmbp <- data.frame("Chr"=c( paste0("chr", as.character(c(1:22))),"chrX","chrY"), "Mbp"=c(249,237,192,183,174,165,153,135,132,132,132,123,108,105,99,84,81,75,69,63,54,57,141,60))}
  tot <- aggregate(list("Total_Reps"=locus$Chr), by = list("Chr"=locus$Chr), FUN = length)
  ave <- aggregate(list("Mean_Units"=(df$All1+df$All2)/2), by = list("Chr"=df$Chr), FUN = mean)
  med <- aggregate(list("Med_Units"=(df$All1+df$All2)/2), by = list("Chr"=df$Chr), FUN = median)
  status <- subset(locus, Status != "STABLE")
  instab <- aggregate(list("Polymorphic_Reps"=status$Status), by = list("Chr"=status$Chr), FUN = length)
  compile <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Chr"), list(chrmbp, tot, ave, med, instab))
  compiple <- compile[is.na(compile)] <- 0
  compile$Mean_Units <- round(compile$Mean_Units)
  compile$Poly_Ratio <- round(compile$Polymorphic_Reps / compile$Total_Reps, 3)
  compile$Reps_per_Mbp <- round(compile$Total_Reps / compile$Mbp, 3)
  compile$Poly_per_Mbp <- round(compile$Polymorphic_Reps / compile$Mbp, 3)
  compile <- compile[order(match(compile$Chr, chrmbp$Chr)),]
  write.table(compile, paste0(output, "_by_chr.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}




