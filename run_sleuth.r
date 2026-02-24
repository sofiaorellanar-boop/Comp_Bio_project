library(sleuth)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
output_filename <- args[1]

#read in the table you made describing samples and kallisto output
#assign to variable name stab 
stab = read.table("metadata.txt",header = TRUE)
stab$path <- file.path("results", stab$sample)

#initialize sleuth object using sleuth_prep function from sleuth library
so = sleuth_prep(stab)

#fit a model comparing the two conditions 
so = sleuth_fit(so, ~condition, 'full')

#fit the reduced model to compare in the likelihood ratio test 
so = sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test for differential expression between conditions 
so = sleuth_lrt(so, 'reduced', 'full')

#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant <- sleuth_table %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::arrange(pval) %>%
  dplyr::select(target_id, test_stat, pval, qval)

#print top 10 transcripts
head(sleuth_significant, n=10)

#write FDR < 0.05 transcripts to file
write.table(sleuth_significant, file = output_filename,quote = FALSE,row.names = FALSE,sep = "\t")
