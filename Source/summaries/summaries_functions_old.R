
by_locus = function(df, output, savename, build){
  df$avealleles <- (df$All1+df$All2)/2
  avera <- aggregate(list("Mean"=df$avealleles), by = list("Call_ID"=df$Call_ID), FUN = mean)
  avera <- avera[order(avera$Call_ID),]
  
  mdata <- melt(df[c("Call_ID", "All2", "All2")], id="Call_ID")
  medi <- aggregate(list("Median" = mdata$value), by=list("Call_ID" = mdata$Call_ID), median)
  medi <- medi[order(medi$Call_ID),]
  
  df <- merge(df,medi,by="Call_ID")
  df$stab <- ifelse(df$All1!=df$Median & df$All2!=df$Median, 1, 
                    ifelse(df$All1==df$Median & df$All2!=df$Median, 0.5, 
                           ifelse(df$All1!=df$Median & df$All2==df$Median, 0.5, 0)))
  
  ref <- aggregate(list("Ref_Units"=df$Ref_Units), by = list("Call_ID"=df$Call_ID, "Chr"=df$Chr, "Start"=df$Start), FUN = mean)
  counts <- aggregate(list("Count"=df$Call_ID), by = list("Call_ID"=df$Call_ID), FUN = length)
  stability <- aggregate(list("Stab"=df$stab), by = list("Call_ID"=df$Call_ID), FUN = mean)
  maxreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("Call_ID"=df$Call_ID), FUN = max)
  minreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("Call_ID"=df$Call_ID), FUN = min)
  stddev <- aggregate(list("Alls"=c(df$All1,df$All2)), by = list("Call_ID"=c(df$Call_ID,df$Call_ID)), FUN=sd)
  
  counts <- counts[order(counts$Call_ID),]
  ref <- ref[order(ref$Call_ID),]
  stability <- stability[order(stability$Call_ID),]
  maxreps <- maxreps[order(maxreps$Call_ID),]
  minreps <- minreps[order(minreps$Call_ID),]
  stddev <- stddev[order(minreps$Call_ID),]
  
  compile <- data.frame("Call_ID"=ref$Call_ID, "Chr"=ref$Chr, "Start"=ref$Start, "Ref_Units"=ref$Ref_Units, "Mean_Units"=round(avera$Mean), "Med_Units"=round(medi$Median), "Min_Rep"=pmin(minreps$All1, minreps$All2), "Max_Rep"=pmax(maxreps$All1, maxreps$All2), "Hits"=counts$Count, "Standard_Dev"=stddev$All,"Instability_Rating"=stability$Stab)
  compile$Status <- ifelse(compile$Instability_Rating == 0, "STABLE", "UNSTABLE")
  
  anno_dir <- paste0(output, "/Annonvar")
  system(paste0("mkdir ", anno_dir))  
  anno_out <- paste0(anno_dir, "/", savename)
  anno <- compile[c("Chr", "Start")]
  anno$End <- (anno$Start)+(3*compile$Ref_Units-1)
  anno$RA <- 0
  anno$AA <- "-"
  write.table(anno, paste0(anno_out, "_Anno.csv"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  build <- paste0("hg", build)
  system(paste0("perl /home/dannear/Binaries/annovar/annotate_variation.pl -out ", anno_dir, "/", savename,  "_Annovar -build ", build ," ", anno_out, "_Anno.csv /home/dannear/Binaries/annovar/humandb"))
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
  write.table(compile, paste0(output, savename, "_Locus_summary.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")  
  return(compile)
}


ratio_of_stability = function(df, output) {
  df <- df[c("Call_ID", "Mean_Units", "Med_Units", "Instability_Rating", "Status")]
  df$Med_Units <- round(df$Med_Units)
  unstable_reps <- subset(df, df$Status != "STABLE")
  stable_reps <- subset(df, df$Status == "STABLE")
  count_stable <- aggregate(list("Stable_Reps"=stable_reps$Med_Units), by = list("Med_Units"=stable_reps$Med_Units), FUN = length)
  count_unstable <- aggregate(list("Unstable_Reps"=unstable_reps$Mean_Units), by = list("Med_Units"=unstable_reps$Med_Units), FUN = length)
  count_total <- merge(count_stable, count_unstable, all = TRUE)
  count_total[is.na(count_total)] <- 0
  count_total$Total_Reps <- count_total$Stable_Reps+count_total$Unstable_Reps
  count_total$Flagged_Unstable <- (count_total$Unstable_Reps)/count_total$Total_Reps
  
  write.table(count_total, paste0(output, "_by_heterozygosity.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}

by_chr = function(df, locus, output){
  chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  mbp <- c(249,237,192,183,174,165,153,135,132,132,132,123,108,105,99,84,81,75,69,63,54,57,141,60) 
  chrmbp = data.frame("Chr"=chr, "Mbp"=mbp)
  tot <- aggregate(list("Reps"=locus$Chr), by = list("Chr"=locus$Chr), FUN = length)
  tot$Chr <- ifelse(tot$Chr=="X", 23, ifelse(tot$Chr=="Y", 24, tot$Chr))
  tot <- tot[order(tot$Chr),]
  ave <- aggregate(list("Ave"=(df$All1+df$All2)/2), by = list("Chr"=df$Chr), FUN = mean)
  ave$Chr <- ifelse(ave$Chr=="X", 23, ifelse(ave$Chr=="Y", 24, ave$Chr))
  ave <- ave[order(ave$Chr),]
  med <- aggregate(list("Med"=(df$All1+df$All2)/2), by = list("Chr"=df$Chr), FUN = median)
  med$Chr <- ifelse(med$Chr=="X", 23, ifelse(med$Chr=="Y", 24, med$Chr))
  med <- med[order(med$Chr),]
  status <- subset(locus, Status != "STABLE")
  instab <- aggregate(list("Status"=status$Status), by = list("Chr"=status$Chr), FUN = length)
  instab <- instab[order(instab$Chr),]
  z <- data.frame("Chromosome"=chrmbp$Chr, "Number_of_Repeats"=tot$Rep, "Mean_Rep_Len"=ave$Ave, "Median_Rep_Len"=med$Med, "Repeats_per_Mbp"=round(tot$Reps/chrmbp$Mbp, 3), "Poly_per_Mbp"=round(instab$Status/chrmbp$Mbp, 3), "Stable_Reps"=tot$Reps-instab$Status, "Unstable_Reps"=instab$Status)
  
  write.table(z, paste0(output, "_by_chromosome.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}


by_sample = function(df, locus, output){

  df <- merge(df, locus[c("Call_ID", "Med_Units", "Mean_Units", "Standard_Dev")], by="Call_ID")
  
  poly <- subset(df, df$All1 != Med_Units | df$All2 != Med_Units)
  pa <- aggregate(poly$Sample_ID, by = list(poly$Sample_ID), FUN = length)
  a <- aggregate(df$Sample_ID, by = list(df$Sample_ID), FUN = length)
  b <- aggregate((All1+All2)/2 ~Sample_ID, data = df, FUN = mean)
  mdata <- melt(df[c("Sample_ID", "All2", "All2")], id="Sample_ID")
  c <- aggregate(mdata$value, by = list(mdata$Sample_ID), FUN = max)
  z <- data.frame("Sample_ID"=a$Group.1, "Highest_Repeat"=c$x, "Total_Repeats"=a$x, "Mean_Rep_length"=round(b$"(All1 + All2)/2", 2), "Unstable_Reps"=pa$x, "Percentage_Unstable_Reps"=round((pa$x/a$x*100), 2))
  z <- z[order(z$Sample_ID),]
  
  write.table(z, paste0(output, "_by_sample.csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}

