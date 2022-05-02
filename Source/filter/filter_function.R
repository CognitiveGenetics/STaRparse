#!/usr/bin/env Rscript

filter_by_coverage <- function(reads, coverage){
	# RENAME COlUMNS
        names(reads) <- c("Call_ID", "Sample_ID", "Chr", "Start", "End", "GT", "Ref_Units", "Allele1_Units", "Allele2_Units")
                
        # MAKE JSON DATA ACCESSIBLE
        coverage$MaxFlankingRead <- str_replace_all(coverage$MaxFlankingRead, "\\(", "")
        coverage$MaxFlankingRead <- str_replace_all(coverage$MaxFlankingRead, "\\)", "")
        coverage$MaxFlankingRead <- str_replace_all(coverage$MaxFlankingRead, " ", "")
        coverage$MaxFlankingRead <- lapply(str_split(coverage$MaxFlankingRead, ','), as.integer)
        coverage$MaxFlankingRead <- lapply(coverage$MaxFlankingRead, max)
                
        # REMOVE NON / NA READS
        reads <- reads[complete.cases(reads), ]
        reads <- reads[!(reads$Allele1_Units=="." | reads$Allele2_Units=="."),]
                
        # MERGE READ AND COVERAGE DATAFRAMES
        df2 <- merge(reads, coverage, x.by=c("Call_ID", "Sample_ID"))
                
        # DELETE NON-MERGED DATAFRAMES
        rm(reads, coverage)
                
        # FILTER OUT BASED ON POOR COVERAGE
        df2 <- subset( df2,  !( df2$MaxSpanningRead == "()" & df2$MaxInrepeatRead == "()" & df2$MaxFlankingRead < df2$Ref_Units ), 
		select =  c("Call_ID", "Sample_ID", "Chr", "Start", "End", "GT", "Ref_Units", "Allele1_Units", "Allele2_Units") )
                
        # RETURN FILTERED DATA FRAME
        return(df2)
        }
