# hapeQTL

`hapeQTL` is a package containing two functions that calculates statistical associations between the level of expression at allele-specific expression (ASE) sites and the presence of *cis* sequence variants by running permutations of the input data and applying a generalized linear model (GLM) each time.


## Dependencies 
* This package was built in R (>= 3.5.0) and requires the following packages to run:
* AER
* biomaRt
* GenomicRanges
* ggplot2
* ggpubr
* plyr
* dplyr
* purrr
* broom
* tidyr
* reshape2
* rtracklayer
* stringr

## Input files 
Only those columns that are required are shown. The sample, legend and hap files are in the IMPUTE format (https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html).

Both the legend and sample files have headers, but the hap does not. Your hap file should have the same number of rows as that of the legend (once the legend's header is accounted for) and should also have twice the number of columns as the sample has rows, because each individual will have two haplotypes. Below are examples of acceptable formats. 

The ASE file must contain the positions of the ASE sites, including the chromosome they are on, the individuals in which they are found, and the alleles at the different haplotypes, including their respective read counts.


### Format of samples file
One column must contain the individual IDs, another, their sex and the third, the population to which they belong.
```
  ID        SEX    POP
 HG00096    male    GBR
 HG00097  female    GBR
 HG00099  female    GBR
 NA20827    male    TSI
```

### Format of legend file
One column must contain the ID of the polymorphisms, another their position, and the final three the reference allele, alternative allele and type of polymorphism involved.
```
 ID                     pos     allele0  allele1
10:60515:C:T            60515     C        T
rs148087467:60523:T:G   60523     T        G
rs147855157:61372:A:C   61372     A        C
```

### Format of haplotypes file
Here, the 0 refers to the reference allele, and the 1, the alternative.
```
0 0 1 0 0 0 0 0 0 0
0 0 0 0 1 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0
```

### Format of ASE file
One column for the chromosome, another for the position, another for the individual ID, two for the reference and alternative alleles, along with two more columns for their corresponding read counts.
```
 chr    start   end     ref     alt     Ind       refCount    altCount
1       135031  135032  G       A       HG00276   19          0
1       135031  135032  G       A       HG00282   12          0
1       135031  135032  G       A       NA11831   10          0
1       135031  135032  G       A       NA19093   12          0

```


## Functions
1. `Gen.input`is a function to generate an input for the model to run. Separating out into two functions saves memory and greatly speeds up parallelisation for what is a computationally demanding task in the second function. The input files are formatted, gene and transcript start site information added to the ASE file, and finally output in .RData format to be loaded into the second function.
2. `Run.Model` is a function that allows you to measure statistical association between nearby regulatory variants and the level of expression at a heterozygous coding polymorphism, controlling for factors such as sex and population, by utilising a generalized linear model and applying permutations to the data in order to provide a robust p-value. As the function supports parallelisation, a number of .txt files equal to the number of tasks, `numTasks`, specified in the function will be outputted.  


## Installing in R ##

### First, load the devtools library: 
```
library(devtools)
```

### Then, install the package from GitHub:
```
install_github("ac1990/hapeQTL")
```

### Then load package as normal:
```
library("hapeQTL")
```

### Download the sample files in R
```
download.file("https://www.dropbox.com/s/cosadgz59hmfoxx/sample_ASE.txt?dl=1","sample_ASE.txt")
download.file("https://www.dropbox.com/s/eycoqh5s6bh8lyq/sample_legend.leg?dl=1","sample_legend.leg")
download.file("https://www.dropbox.com/s/e8sfx3gbziserkh/sample_haplotypes.hap?dl=1","sample_haplotypes.hap")
download.file("https://www.dropbox.com/s/8su6scg0ojkc8pm/sample_Samples.txt?dl=1","sample_Samples.txt")
```

### Or download manually via a web browser (make sure they're in your working directory, default set to PackageTestWork)
* https://www.dropbox.com/s/cosadgz59hmfoxx/sample_ASE.txt?dl=1
* https://www.dropbox.com/s/eycoqh5s6bh8lyq/sample_legend.leg?dl=1
* https://www.dropbox.com/s/e8sfx3gbziserkh/sample_haplotypes.hap?dl=1
* https://www.dropbox.com/s/8su6scg0ojkc8pm/sample_Samples.txt?dl=1


## `Gen.input`
This function allows you to output .RData files to be used as input for the analysis run in the second function. This only needs to be done once for each chromosome which ensures the running of the model is not too computationally demanding. The processing of input files need not be done repeatedly, as they are constants for whichever parameters you wish to set in `Run.Model.R`. Hence separating this processing out allows tweaks to be made to the parameters in the second function, without having to re-run the processing each time.
       
### Arguments:

1. ASE_file: Specify file containing allele-specific expression information.
2. legend_file: Specify legend file containing SNP IDs and position info, etc.
3. haplotypes_file: Provide haplotypes corresponding to SNPs in legend file, for given individuals found in samples file.
4. samples_file: File containing cohort information.
5. output_prefix: Pre-fix to the name of the output file.
6. species: Select same species as cohort. Defaults to "hsapiens".
7. ensembl_version: Specify version of Ensembl to download from. Defaults to NULL, which is the most recent build.

### Examples:
#### Annotating the input files using the most recent human gene set
```
Gen.input('path_to_ASE_file/ASEfile.txt.gz',
                     'path_to_legend_file/file.legendfile.gz',
                    'path_to_haplotypes_file/hapfile.hap.gz',
                    'path_to_samples_file/samples.txt')
```                     
#### Annotating the input files using a previous human gene set (the same as the one that corresponds to the sample data
```
Gen.input('path_to_ASE_file/ASEfile.txt.gz',
                     'path_to_legend_file/file.legendfile.gz',
                    'path_to_haplotypes_file/hapfile.hap.gz',
                    'path_to_samples_file/samples.txt', ensembl_version=78)
```     
#### Annotating the input files using the most recent mouse gene set
```    
Gen.input('path_to_ASE_file/ASEfile.txt.gz',
                     'path_to_legend_file/file.legendfile.gz',
                   'path_to_haplotypes_file/hapfile.hap.gz',
                     'path_to_samples_file/samples.txt',
                     '                     'mmusculus')
```

The RData files outputted from the function, Gen.input, can now be used to run the analysis in the second function, below.


## `Run.model`            
### Arguments: 
1. inputObj: This is the RData file outputted from the first function, Gen.input
2. output_prefix: This is the pre-fix of the name of the output file
3. task: This analysis is very computationally burdensome. To speed up the process it is an advantage to split it up into tasks that may be run on multiple nodes, concurrently. 
4. totalTasks: Set the total number of tasks to split the analysis into. Defaults to 100.
5. minInd: Specify the minimum number of individuals in which an ASE site must be found before in order to be included in the analysis. Defaults to 10
6. numPerms: Select how many permutations to run. Along with splitting the process up into simultaneous tasks, this is the biggest factor in determining how long the analysis will take. However, the more permutations, in general and up until a point, the more precise and accurate the results may be; for example, if set to 100, the minimum p-value that can possibly be reached as a result of permutations, is 0.01. Defaults to 100,000. 
7. TSSwindow: The transcript start site window. This represents the distance over which nearby variants will be selected,  either side of the transcript start site. Defaults to 500kb 
8. pval_threshold: There is a theoretical minimum p-value for each particular combination of reference and alternative alleles for a given set of individuals for a given nearby variant of an ASE site. This parameter sets the upper limit. Default is 0.00005. 
9. other_all: Specify whether or not to account for the allele on the haplotype to which the expressed site does not belong. Defaults to FALSE
10. seed: Specify the number that the seed should be set to. The seed is the starting point used in the generation of a sequence of random numbers. Defaults to 10
        
          
### Examples:
#### Run model with task set to 10, for 100,000 permutations, a transcript start site window of 500kb and a theoretical p-value threshold of 0.00005
```
Run.Model('input_file.RData', 10, 
        'output_prefix')
```
#### Run model with task set to 10, for 10,000 permutations, a transcript start site window of 500kb and a theoretical p-value of 0.00005
```
Run.Model('input_file.RData', 10,
         'output_prefix',
         numPerms=10000, pval_threshold=10000)
 ```    
#### Run model with task set to 2, for 10,000 permutations, a transcript start site window of 1Mb and a theoretical p-value of 0.00005
```
Run model with task\tset to 2, for 10,000 permutations, a transcript start site window of 1Mb and a theoretical p-value of 0.00005
#' Run.Model('input_file.RData', 2,
#'         'output_prefix',
#'         numPerms=10000,TSSwindow=100000, pval_threshold=10000)
 ```    

NB: The smallest possible p-value attainable as a result of running permutations is 1/numPerms. Hence, there is no advantage to setting the minimum p-value threshold to below this number.


## Output format

The first four columns are the ID of the ASE site, its position, the individual, and the ID of the nearby variant. Then follows 7 columns representing the test statistics of the intercept, reads, sex, sub-populations and variant, with the next 7 columns representing their corresponding p-values. The next two columns are the transcript start site of the ASE site and the gene to which the nearby variant belongs. Where the nearby variant is not found in a gene, the ID of the variant itself is given. The final two columns are the number of permutations, and the number of permutations in which the test statistic is exceeded.  





## Selecting significant sites
To calculate the permuted p-values, in R do the following to your results dataframe, `df`:
```
#### Divide numPermExceed column by the numPerm column
df2 <- df[which(as.numeric(df$numPerm) > 0),]
df2$permuted_p <-as.numeric(df2$numPermExceed)/as.numeric(df2$numPerm)
```

## Plotting functions
You can use the `plotASEMetrics` function to visualise some of your ASE data. It shows how many individuals are heterozygote at each ASE site in the first plot. In the second plot, it shows how many genes contain varying number of heterozygous sites. The third plot provides a survey of how many sites show a certain proportion of reads that carry the reference allele rather than the alternative allele, which is effectively a check for any potential allelic mapping bias. The fourth plot shows the extent of ASE across the individual sites and individuals, plotting allele ratio against the total number of reads and colouring the points according to their binomial p-value. Simply load your output from the Gen.input function as follows:

`plotASEMetrics(aseDat)`

You can also use the `plotQTL` function to construct violin plots which show the distribution of the expression values on the basis of the status of the nearby allele, on both haplotypes. It requires three inputs be specified inside the function. First, the output from the Gen.input function; secondly, an ASE site; thirdly, a putative cis-regulatory variant, which you can find in the output of the Run.model function. An example is given below:

`plotQTL(aseDat, "rs9306160", "rs2838351", otherAll = TRUE)`



### Having difficulty installing some packages?
```
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
biocLite("GenomicRanges")
```


EOF
