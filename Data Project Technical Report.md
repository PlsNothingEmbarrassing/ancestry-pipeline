# Technical Report: Pipeline to Analyse both Ancestry Similarity and Sex from Genotype Data

## Background

Genetic epidemiology is the study of the relationship between genetic variation and the . As our capacity to store and process data has grown, increasingly rigorous quality control (QC) measures are necessary to minimise bias and error in results. Two aspects of genomic QC that I will be focussing on are ancestry inference and sex determination. 
### Ancestry Inference
Ancestry inference is a key component in QC as population structure can affect genetic association studies. Failing to account for differences in genetic ancestry in cases and controls in genome wide association studies (GWAS) can lead to false positive results (Byun J et al,. 2017). Self reported ancestry can often be inaccurate due a number of factors such as lack of knowledge, differences in perceived ancestry and reality, and social and cultural influences. One common method of inferring ancestry is through Principal Component Analysis (PCA), which is a dimensionality reduction method which is used to produce several uncorrelated variables (or PCs)  from a data matrix containing observations across a number of potentially correlated variables (Anderson, CA et al., 2010). By projecting individual genotypes onto PCs, samples can be clustered based on genetic similarities.
Using labelled reference datasets can improve the reliability of the ancestry estimation due to being able to place samples within pre-established clusters.

### Sex Determination
Sex determination is also an important QC step as discrepancies between reported and genetic sex can indicate issues with the data, such as plating errors and sample mix ups (Anderson, CA et al., 2010). Genetic sex is typically derived from evaluating heterozygosity on the X chromosome as is the case with the PLINK "check-sex" method (www.cog-genomics.org, n.d.). Checking X chromosome heterozygosity can also reveal chromosome abnormalities  such as Kleinfelter syndrome (male XXY), Turner syndrome (females X0), or mosaic individuals (Turner, S et al., 2011).
## Methods

## Data
### 1000 Genomes Project Data
The 1000 genomes project data is used as the reference dataset for my ancestry pipeline. The data is comprised of 2504 individuals from 26 populations representing a range of global populations.
#### Data Components
1. **Population files:** Population file contains individual ID (IID), family ID (FID), and Population label e.g. GBR, IBS.
2. **Superpopulation file:** Superpopulation file contains individual ID (IID), family ID (FID), and Superpopulation label e.g. EUR, SAS.
3. **Genotype Data:** Genotype data is stored in .bed, .bim and .fam files. Reported sex is stored in .fam file. Data contains SNP data along with other variant types.
4. **Ancestry Informative Markers (AIMs):** AIMs data contains list of SNPs that have been selected as Ancestry Informative Markers. This specifies which SNPs the data should be limited to so that the 
## Software
The code for this project was split into 4 sections:
1.  Defining ancestry metrics and performing PCA to extract PC data using PLINK.
2. Evaluating results of PCA and identifying individuals in the test dataset who show similar ancestry to GBR reference individuals in R.
3. Using PLINK to identify a subset of individuals in the test dataset whose genotypes identify them as female and comparing to the reported sex values.
4. Evaluating the results of the sex-check and generating appropriate visualisations in R.
### Defining Ancestry Metrics
#### Soft Coding and Validation Methods
For the first part of the task, I developed a bash script in which the user could specify the paths for the reference data directory, test data directory, path to linkage disequilibrium (LD) ranges, output and results. I did this using optargs which allows the user to pass arguments to the script when invoking through the CLI. 

```bash
# Default values for ease of use for myself.
TEST=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results
REFERENCE=/home/c.c24082291/PLINK_Bioinfo/database/genotypes/all-1000g-phase3-chrall-mac5-v2/all-1000g-phase3-chrall-mac5-v2
HIGHLD=/home/c.c24082291/PLINK_Bioinfo/database/ranges/high-ld-b37.ranges
OUTDIR=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/tmp
RESULTS=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results

# Parse command line arguments
while getopts ":t:r:h:o:s:" opt; do

    case ${opt} in

        t ) TEST=$OPTARG ;;

        r ) REFERENCE=$OPTARG ;;

        h ) HIGHLD=$OPTARG ;;

        o ) OUTDIR=$OPTARG ;;

        s ) RESULTS=$OPTARG ;;

        * ) usage ;;

    esac

done
```
*Code block 1*

To avoid confusion of how to pass arguments properly, I added a "usage" function that would run if incorrect flags were passed on execution. This would output a help message which would inform a user how to structure their command. 
```bash
# Function to display help message

usage() {

    echo "Usage: $0 -t TEST -r REFERENCE -h HIGHLD -o OUTDIR -s RESULTS"

    echo "  -t TEST       Path to the test dataset"

    echo "  -r REFERENCE  Path to the reference dataset"

    echo "  -h HIGHLD     Path to the high LD ranges file"

    echo "  -o OUTDIR     Path to the output directory"

    echo "  -s RESULTS    Path to the results directory"

    exit 1

}
```
*Code block 2*

Additionally I added simple data validation steps to ensure that the user provided arguments for each necessary flag.
```bash
# Check if all required arguments are provided

if [ -z "$TEST" ] || [ -z "$REFERENCE" ] || [ -z "$HIGHLD" ] || [ -z "$OUTDIR" ] || [ -z "$RESULTS" ]; then

    usage

fi
```
*Code block 3*
#### Main Script body

##### Extracting SNP IDs from both datasets

Now that the files necessary files and directories have been set, we want to start by identifying the SNP IDs in both the test and reference dataset. In the data used in this task that is stored in the 2nd column in the .bim file from both datasets, so using awk we can extract these and store separately. I also print the number of SNPs found in each to the console as an easy way to see if there are errors in reading the files e.g. 0 SNPs detected in either sample.

```bash
module load plink/1.9

# Get SNPs from both datasets

echo "awking snps"

awk '{print $2}' "${REFERENCE}.bim" | sort > "${OUTDIR}/reference.snps"

awk '{print $2}' "${TEST}.bim" | sort > "${OUTDIR}/test.snps"

  
# Output number of snps in each dataset
echo "Number of snps in reference dataset: $(wc -l < "${OUTDIR}/reference.snps")"

echo "Number of snps in test dataset: $(wc -l < "${OUTDIR}/test.snps")"
```
*Code block 4*

##### Finding SNP IDs Common to both Datasets

 In the next step, we use the `comm` command with the `-12` flag to compare the two files and keep only lines which are common to both test and reference data. These common SNP IDs are stored in an `.extract` file. The number of common SNP IDs are printed to the console to inform the user of any unexpected absences of output. I could have added further validation to halt script execution if there were not any common SNPs.
 
 ```bash
 # Use comm to identify common snps

# Use -12 flag to suppress lines unique to each file

comm -12  \

"${OUTDIR}/reference.snps" \

"${OUTDIR}/test.snps" \

> "${OUTDIR}/common_to_both.extract"

  

echo "Number of snps in common: $(wc -l < "${OUTDIR}/common_to_both.extract")"
```
*Code block 5*

##### Extracting SNP data matching the shared SNP IDs

The next step is to extract the data rows with the SNP information that are common to both reference and test datasets. This is different to the previous steps as we are extracting all of the related information not just the SNP ID. The common SNP data from each dataset is stored in a new `.bim` file for further processing.

```bash
# limit to commons snps

# Use plink to extract common snps from ref dataset

# Use --extract flag to specify snps to keep

plink \

--bfile "${REFERENCE}" \

--extract "${OUTDIR}/common_to_both.extract" \

--make-bed \

--out "${OUTDIR}/reference_common"
```
*Code block 6*

##### Remove Variants with low Minor Allele Frequency (MAF)

The next step in the pipeline is to remove any SNPs with very low MAF. This is important as association signals from these SNPs are less informative as they occur in only a small number of individuals (Anderson, CA et al., 2010). The filtered dataset is then stored for further processing.

```bash
# Use plink to extract common snps from in the new ref dataset

# Use maf flag to specify minor allele frequency

plink \

--bfile "${OUTDIR}/reference_common" \

--maf 0.05 \

--make-bed \

--out "${OUTDIR}/reference_common_maf5"

echo "Number of snps in common with MAF > 0.05: $(wc -l < "${OUTDIR}/reference_common_maf5.bim")"
```
*Code block 7*

##### Remove Ambiguous SNPs

The next process is to remove ambiguous SNPs. This is where the alleles are complementary base pairs e.g. A/T, C/G.
These are considered ambiguous because without knowing the strand the same SNP could be reported as forward or reverse. For example a SNP reported as A/T on the forward strand could be T/A on the reverse strand. Without knowing the strand the true allele configuration is ambiguous and can lead to errors so it is important to exclude them from analysis.

```bash
# Remove ambiguous snps

# Check if the snp is ambiguous

# Use cut to extract the snp id and write to a file

awk \

'($5=="T"&&$6=="A")||($5=="A"&&$6=="T")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C")' \

"${OUTDIR}/reference_common_maf5.bim" | cut -f 2 > "${OUTDIR}/ambiguous.exclude"

# Use plink to exclude ambiguous snps from the ref dataset

plink \

--bfile "${OUTDIR}/reference_common_maf5" \

--exclude "${OUTDIR}/ambiguous.exclude" \

--make-bed \

--out "${OUTDIR}/reference_common_maf5_noWS"
```
*Code block 8*



##### Limit to Ancestry Informative Markers (AIMs)

The next step is to limit the filtered SNPs to ancestry informative markers which are SNPs that have large frequency divergences among populations from different geographic regions and are can be used to estimate an individuals ancestry (Felkl, A.B et al., 2023). This also simplifies further analysis by focussing on a smaller set of markers which have the most importance.

```bash
# Limit to ancestry informative markers

# Use plink to extract ancestry informative markers from the ref dataset

plink \

--bfile "${OUTDIR}/reference_common_maf5_noWS" \

--extract "${REFERENCE}.aims" \

--make-bed \

--out "${OUTDIR}/reference_common_maf5_noWS_aims"
```
*Code block 9*


##### Identify and Remove Regions of High Linkage Disequilibrium

We now need to use plink to identify SNPs which overlap the ranges which are flagged as high LD. LD is a statistical association between variants meaning that SNPs are strongly correlated with each other which can lead to bias and, if presuming SNP independence, could cause errors. Pruning these regions ensures that the analysis is not dominated by a small number of highly correlated samples.

Using plink we can create a set of SNPs that overlap with the ranges specified in our `ranges` file which are high LD regions. The SNPs are stored so they can be excluded in further analysis.

```bash
# Remove regions of high LD

# Use plink to identify set of snps that overlap the ranges flagged as high LD

plink \

--bfile "${OUTDIR}/reference_common_maf5_noWS_aims" \

--make-set "${HIGHLD}" \

--write-set \

--out "${OUTDIR}/reference_common_maf5_noWS_aims"

  

# Use plink to exclude snps in the high LD regions

plink \

--bfile "${OUTDIR}/reference_common_maf5_noWS_aims" \

--exclude "${OUTDIR}/reference_common_maf5_noWS_aims.set" \

--make-bed \

--out "${OUTDIR}/reference_common_maf5_noWS_aims_noLD"
```
*Code block 10*

Next we use plink to generate a list of SNPs in a `.prune.in` file to prune from our datasets.

```bash
# LD prune

# Use plink to generate a LD pruned list of snps .prune.in

plink \

--bfile "${OUTDIR}/reference_common_maf5_noWS_aims_noLD" \

--indep-pairwise 50 10 0.1 \

--out "${OUTDIR}/ld_independent"

# Specify the parameters for LD pruning (window size, step size, r^2 threshold)
```
*Code block 11*

##### Extract LD Limited SNPs from Reference and Test Dataset

Using plink we can now extract the SNP information from the LD limited SNPs in the test and reference dataset in preparation to merge both.

```bash
# LD limit and merge

# Use plink to extract the LD limited SNPs from the ref and test dataset

plink \

--bfile "${REFERENCE}" \

--extract "${OUTDIR}/ld_independent.prune.in" \

--make-bed \

--out "${OUTDIR}/reference_premerge"

  

# Apply to the test dataset

plink \

--bfile "${TEST}" \

--extract "${OUTDIR}/ld_independent.prune.in" \

--make-bed \

--out "${OUTDIR}/test_premerge"
```
*Code block 12*

##### Merge Datasets

The next step is to merge the two processed datasets. We begin by doing an initial merge attempt. If there are any strand alignment issues, the merge will fail and plink generates a list of SNPs which are problematic in a `.missnp` file.

```bash
# Merge datasets

plink \

--bfile "${OUTDIR}/reference_premerge" \

--bmerge "${OUTDIR}/test_premerge" \

--make-bed \

--out "${OUTDIR}/combined"
```
*Code block 13*

If the merge fails, the next step is to flip the strand of the SNPs with strand alignment issues and attempt the merge again.

```bash

plink \

--bfile "${OUTDIR}/reference_premerge" \

--flip "${OUTDIR}/combined-merge.missnp" \

--make-bed \

--out "${OUTDIR}/reference_flipped"

  

# Use plink to merge the datasets including those with strand alignment issues

plink \

--bfile "${OUTDIR}/reference_flipped" \

--bmerge "${OUTDIR}/test_premerge" \

--make-bed \

--out "${OUTDIR}/combined"
```
*Code block 14*

##### Perform PCA on merged and processed dataset

The penultimate step is to calculate principal components from the merged dataset. The output gives us an `.eigenvec` and `.eigenval` file with values from the matrix computation from the PCA.

```bash
# Calculate principal components

# Use plink to calculate the principal components

plink \

--bfile "${OUTDIR}/combined" \

--pca \

--out "${RESULTS}/combined"
```
*Code block 15*

##### Clean Up Temporary Files

The final step is to remove intermediatory files generated throughout the pipeline. This of course can be left out or modified to retain certain files if necessary. 

```bash
# Clean up temporary files

rm "${OUTDIR}/reference.snps" "${OUTDIR}/test.snps" "${OUTDIR}/common_to_both.extract"

rm "${OUTDIR}/reference_common.bed" "${OUTDIR}/reference_common.bim" "${OUTDIR}/reference_common.fam"

rm "${OUTDIR}/reference_common_maf5.bed" "${OUTDIR}/reference_common_maf5.bim" "${OUTDIR}/reference_common_maf5.fam"

rm "${OUTDIR}/ambiguous.exclude"

rm "${OUTDIR}/reference_common_maf5_noWS.bed" "${OUTDIR}/reference_common_maf5_noWS.bim" "${OUTDIR}/reference_common_maf5_noWS.fam"

rm "${OUTDIR}/reference_common_maf5_noWS_aims.bed" "${OUTDIR}/reference_common_maf5_noWS_aims.bim" "${OUTDIR}/reference_common_maf5_noWS_aims.fam"

rm "${OUTDIR}/reference_common_maf5_noWS_aims.set"

rm "${OUTDIR}/reference_common_maf5_noWS_aims_noLD.bed" "${OUTDIR}/reference_common_maf5_noWS_aims_noLD.bim" "${OUTDIR}/reference_common_maf5_noWS_aims_noLD.fam"

rm "${OUTDIR}/ld_independent.prune.in" "${OUTDIR}/ld_independent.prune.out" "${OUTDIR}/ld_independent.prune.log"

rm "${OUTDIR}/reference_premerge.bed" "${OUTDIR}/reference_premerge.bim" "${OUTDIR}/reference_premerge.fam"

rm "${OUTDIR}/test_premerge.bed" "${OUTDIR}/test_premerge.bim" "${OUTDIR}/test_premerge.fam"

rm "${OUTDIR}/reference_flipped.bed" "${OUTDIR}/reference_flipped.bim" "${OUTDIR}/reference_flipped.fam"

rm "${OUTDIR}/combined-merge.missnp"
```
*Code block 16*

### Evaluating Results of PCA and Identifying GBR-like Individuals

After getting the results of our PCA in the form of our `.eigenvec` and `.eigenval` files, we next imported the data into RStudio in order to perform some exploratory data analysis and visualisation. 
When that is complete we can move on to defining cluster thresholds to define individuals with GBR-like ancestry in the test dataset.
#### Loading and Formatting PCA Data

The first step after loading our libraries is to import our data into RStudio as a tibble using the `readr` library. As our data is headerless we need to define column names. Some columns need to have their datatype explicitly coerced such as FID and IID due to the front section of the samples having integer looking IDs which cause a type error when reading non-int-like rows included later in the data.
```r
eigenvec <- read_table(
  "E:/Bioinformatics/PCA_ASSESSMENT/combined.eigenvec",
  col_names = c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), # Specify the exact names
  col_types = cols(FID = col_character(), IDD = col_character()) # Use specific column types if needed
)
# Didn't work when reading in so mutating to ensure char dtype on IID col
eigenvec <- eigenvec %>%
  mutate(IID = as.character(IID))
  
# single col data so dont need read_table
eigenval <- read_csv(
  "E:/Bioinformatics/PCA_ASSESSMENT/combined.eigenval",
  col_names = c("Eigenvalue") # Assign the colname"
)
```
*Code block 17*
# 1500ish words

#### Generating Scree Plot from `eigenvals`

After importing our data, we need to calculate explained variance for each component and plot them to create a Scree plot so that we can decide upon the number of principal components to retain. 

```r
# calculate variance + cum variance
var_explained <- eigenval$Eigenvalue / sum(eigenval$Eigenvalue) * 100
cum_var <- cumsum(var_explained)

# Create tibble for plotting
scree_data <- tibble(
  PC = 1:length(eigenval$Eigenvalue),
  VarianceExplained = var_explained,
  CumulativeVariance = cum_var
)
```
*Code block 18*

Next I used `ggplot2` and `plotly` libraries to plot the variance explained for each principal component as well as the cumulative variance. I used `plotly` as it provides 'Show closest data on hover' which I find is helpful to explore values on the plot.

```r
p <- ggplot(scree_data, aes(x = PC)) +
  geom_point(aes(y = VarianceExplained), size = 3, color = "blue") +
  geom_line(aes(y = VarianceExplained), color = "blue", linetype = "solid") +
  geom_point(aes(y = CumulativeVariance), size = 3, color = "red") +
  geom_line(aes(y = CumulativeVariance), color = "red", linetype = "dashed") +
  geom_hline(yintercept = 80, linetype = "dotted", color = "gray") +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_minimal()

int_plot <- ggplotly(p)
int_plot
```
*Code block 19*


## Description of Workflow

## Results and Discussion

## References

- Byun, J., Han, Y., Gorlov, I.P., Busam, J.A., Seldin, M.F. and Amos, C.I., 2017. Ancestry inference using principal component analysis and spatial analysis: a distance-based analysis to account for population substructure. _BMC genomics_, _18_, pp.1-12.
-  Anderson, C.A., Pettersson, F.H., Clarke, G.M., Cardon, L.R., Morris, A.P. and Zondervan, K.T., 2010. Data quality control in genetic case-control association studies. _Nature protocols_, _5_(9), pp.1564-1573.
-  Turner, S., Armstrong, L.L., Bradford, Y., Carlson, C.S., Crawford, D.C., Crenshaw, A.T., de Andrade, M., Doheny, K.F., Haines, J.L., Hayes, G. and Jarvik, G., 2011. Quality control procedures for genome‐wide association studies. _Current protocols in human genetics_, _68_(1), pp.1-19.
-  www.cog-genomics.org. (n.d.). _Basic statistics - PLINK 1.9_. [online] Available at: https://www.cog-genomics.org/plink/1.9/basic_stats.
- Felkl, A.B., Avila, E., Gastaldo, A.Z., Lindholz, C.G., Dorn, M. and Alho, C.S. (2023). Ancestry resolution of South Brazilians by forensic 165 ancestry-informative SNPs panel. _Forensic Science International: Genetics_, [online] 64, p.102838. doi:https://doi.org/10.1016/j.fsigen.2023.102838.

‌
## Appendix

### Acronyms 
QC: Quality Control
LD: Linkage Disequilibrium
GWAS: Genome Wide Association Study
PCA: Principal Component Analysis
PC: Principal Component
AIMs: Ancestry Informative Markers
CLI: Command Line Interface
MAF: Minor Allele Frequency
SNP: Single Nucleotide Polymorphism