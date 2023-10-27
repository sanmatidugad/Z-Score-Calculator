rm(list=ls())    #clear the environment in R. Remove all objects (variables, functions, etc.) from the current workspace.
Sys.setenv(TZ='EDT')

ReQ_Modules = c("tidyr", "data.table", "matrixStats")

for (each in 1:length(ReQ_Modules)) {
  if(ReQ_Modules[each] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_Modules[each])
  }
}

#invisible(lapply(ReQ_Modules, library, character.only = TRUE))

suppressWarnings(library(tidyr))
suppressWarnings(library(data.table))
suppressWarnings(library(matrixStats))

args=commandArgs(trailingOnly = T)

normalized_data = NULL    #normalized data
states_data = NULL    #gene states
gene_counts = NULL    #gene counts file from subsampled output
actionable_data = NULL    #actionable gene list
normal_count = NULL    #number of normal samples
tumor_sample = NULL    #tumor sample name


if ("-i" %in% args ){ #&& "-s" %in% args && "-gc" %in% args && "-n" %in% args && "-t" %in% args && "-a" %in% args ) {    
  normalized_data_index = which(args == "-i")
  states_data_index = which(args == "-s")    #Getting index of the "--input" option
  gene_counts_index = which(args == "-gc")
  normal_count_index = which(args == "-n")
  tumor_sample_index = which(args == "-t")
  actionable_data_index = which(args == "-a")
  
  normalized_data = args[normalized_data_index + 1] 
  states_data =  args[states_data_index + 1] 
  gene_counts = args[gene_counts_index + 1]
  normal_count = as.numeric(args[normal_count_index + 1])
  tumor_sample = args[tumor_sample_index + 1]
  actionable_data = args[actionable_data_index + 1]
}

normalized_data = "test1_normalized.csv"
setwd("/home/data/Dropbox (Genomic Expression)/GEx Bioinformatics/12. OneRNA Secondary Pipeline/Sanmati Normalization Experiments/Z_Analysis/")

#print(paste("Normalized File Name: ", normalized_data))
normalized_input = read.csv(normalized_data, row.names = 1, sep = ",")
#print(paste("Read Normalized File"))
states_input = read.csv(states_data, row.names = 1)
#print(paste("Read Gene States File"))
#print(paste(tumor_sample))

counts_input = read.csv(gene_counts, row.names = 1)
#print(paste("Read Gene Counts File"))

actionable_input = read.csv(actionable_data)
#print(paste("Read Actionable File"))


#head(states_input)
#colnames(normalized_input)
normal_count = 23
dim(normalized_input)

tumor_sample = colnames(normalized_input)[normal_count+1:ncol(normalized_input)]
normalized_input[, c(tumor_sample)] = round(normalized_input[, c(tumor_sample)], 3)    # round up tumor data point to 
normalized_input$Normalized_Per_Million = round((normalized_input[,c(tumor_sample)] * 1000000) / counts_input[c(tumor_sample), c("new_totals")] , 3)
normalized_input$Normal_Mean = round(rowMeans(normalized_input[,1:normal_count]), 3)  # calculating Normal Means
normalized_input$Normal_Median = round(apply(normalized_input[,1:normal_count], MARGIN = 1, median), 4)    #calculating Normal Medians
normalized_input$Normal_SD = round(apply(normalized_input[,1:normal_count], 1, sd) , 3)    #calculating Standard Deviation 
normalized_input$Normal_SD[normalized_input$normal_SD == 0] = 1    # replacing standard deviation of zero by 1
normalized_input$Z_score = round((normalized_input[, c(tumor_sample)] - normalized_input$Normal_Mean) / normalized_input$Normal_SD , 3)   #calculating Z-scores
normalized_input$Foldchange = round(normalized_input[, c(tumor_sample)] / normalized_input$Normal_Mean, 3)    # calculate foldchange values
normalized_input[is.na(normalized_input)] = 0    # replace all NA's and NAN's by 0
normalized_input[normalized_input == Inf] = 0    # replacing Inf by 0
normalized_input[normalized_input == -Inf] = 0

states_tumor = states_input[, c("SYMBOL", tumor_sample), drop = FALSE]    #selecting SYMBOL column and tumor sample column in gene state data
colnames(states_tumor) = c("SYMBOL", "Gene_State")    #renaming tumor sample in gene state data    
normalized_input.1 = merge(normalized_input, states_tumor, by = "row.names", all.x = TRUE)    #merging normalized data table with gene state tumor data
head(normalized_input.1)

#Adding False Positive and False Negative Check
normalized_input.1$Potential_FP_or_FN = NULL
#FALSE POSITIVE: median normal < 10, normalized expression < 10, Aberrantly Called (Very High, High, Low, Very Low)
normalized_input.1[which(normalized_input.1[,c(tumor_sample)] < 10 & normalized_input.1[,c("Normal_Median")] < 10 & normalized_input.1[,c("Gene_State")] != "NORMAL"), c("Potential_FP_or_FN") ] = "+"
#FALSE NEGATIVE: median normal < 10, normalized expression < 10, Not Aberrantly Called (Normal)
normalized_input.1[which(normalized_input.1[,c(tumor_sample)] < 10 & normalized_input.1[,c("Normal_Median")] < 10 & normalized_input.1[,c("Gene_State")] == "NORMAL"), c("Potential_FP_or_FN")] = "-"
normalized_input.1[is.na(normalized_input.1)] = 0

rownames(normalized_input.1) = normalized_input.1[,1]    # changing back rownames to original

## Selecting Columns for Output File
normalized_input.2 = normalized_input.1[, c("SYMBOL", tumor_sample, "Normalized_Per_Million","Normal_Mean", "Normal_Median", "Normal_SD", "Z_score", "Foldchange", "Gene_State", "Potential_FP_or_FN")]    #selecting columns for final output
write.csv(normalized_input.2, paste(tumor_sample, "-z_stats.csv", sep = ""))  

## Only Actionable Genes.
a = actionable_input[trimws(actionable_input$gene) %in% trimws(normalized_input.2$SYMBOL),]
normalized_input.3 = normalized_input.2[normalized_input.2$SYMBOL %in% a$gene, ]
write.csv(normalized_input.3, paste(tumor_sample, "-z_stats-actionable.csv", sep = ""))

