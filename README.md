The following scripts are methods to compute a p-value for genes under the null model of a binomial distribution in a single-cell RNA sequencing dataset as described in the methods section of the manuscript. 

These scripts are written in python and R, with additional relevant details described below.

**Formatting input dataset:**
  The input single-cell RNA sequencing dataset should be formatted as a csv, with columns representing distinct genes and rows representing distinct cells. Thus a value X at row Y and column Z would indicate that X UMI’s are found in cell Y of the transcript counts of gene Z.

**Format of Output:**
  The output of this program will return a csv with two columns. Each row will list a gene in the first column, and said gene’s corresponding p-value in the second column. 

**Python:**

        Required packages:
                Pandas, numpy, scipy

        For ease of execution, the script should be in the same location as the input file, and the 
        variable “large_root” should have the string of the name of the input dataset. 

        The output file will be named using the name of the input dataset, with the phrase “_PVAL.csv” appended to it.


**R**:

        For ease of execution, the script should be in the same location as the input file, and the 
        variable “data” should call the read.csv function with the string of the name of the input dataset. 


The output file name can be selected through the input to the variable output_file.
