data <- read.csv("100_rep_1zheng9_resamp.csv")
output_file <- "100_rep_1zheng9_resamp_RPVALS.txt"
gene_names <- colnames(data)
number_of_cells <- dim(data)[1]
capture_probability <- 0.05

for (n in 2:length(gene_names)){
  gene <- gene_names[n]
  numerical_array <- as.numeric(unlist(data[gene]))
  if (sum(numerical_array) > 0){
    cells_with_mRNAs <-length(which(numerical_array != 0))
    avg_expression <- round(mean(numerical_array), digits = 6)
    mRNA_amount <- avg_expression/capture_probability
    no_mRNA_in_cell <- exp(mRNA_amount * log(1 - capture_probability))
    cells_with_no_mRNAs <- number_of_cells - cells_with_mRNAs
    x <- seq(cells_with_no_mRNAs,number_of_cells - 1,by = 1)
    pval <- round(sum(dbinom(x, size=number_of_cells, prob=no_mRNA_in_cell)), digits = 6)
    final_output <- paste(gene,",",pval)
    write(final_output, file = output_file, append=TRUE)
  }
}

