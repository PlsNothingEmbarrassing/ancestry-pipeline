# Technical Report: Pipeline to Analyse both Ancestry Similarity and Sex from Genotype Data
SN: 24082291
Name: Jude Hutton
## Background

As our capacity to store and process data has grown, increasingly rigorous quality control (QC) measures are necessary to minimise bias and error in results. Two aspects of genomic QC which are key in genome wide association studies that we were tasked to utilise for this assignment are ancestry inference and sex determination. 
### Ancestry Inference
Ancestry inference is a key component in QC as population structure can affect genetic association studies. Failing to account for differences in genetic ancestry in cases and controls in genome wide association studies (GWAS) can lead to false positive results (Byun J et al,. 2017). Self reported ancestry can often be inaccurate due a number of factors such as lack of knowledge, differences in perceived ancestry and reality, and social and cultural influences. One common method of inferring ancestry is through Principal Component Analysis (PCA), which is a dimensionality reduction method which is used to produce several uncorrelated variables (or PCs)  from a data matrix containing observations across a number of potentially correlated variables (Anderson, CA et al., 2010). By projecting individual genotypes onto PCs, samples can be clustered based on genetic similarities.
Using labelled reference datasets can improve the reliability of the ancestry estimation due to being able to place samples within pre-established clusters.

### Sex Determination
Sex determination is also an important QC step as discrepancies between reported and genetic sex can indicate issues with the data, such as plating errors and sample mix ups (Anderson, CA et al., 2010). Genetic sex is typically derived from evaluating heterozygosity on the X chromosome as is the case with the PLINK "check-sex" method (www.cog-genomics.org, n.d.). Checking X chromosome heterozygosity can also reveal chromosome abnormalities  such as Kleinfelter syndrome (male XXY), Turner syndrome (females X0), or mosaic individuals (Turner, S et al., 2011).
## Methods
The code for this project was split into 4 sections:
1.  Defining ancestry metrics and performing PCA to extract PC data using PLINK.
2. Evaluating results of PCA and identifying individuals in the test dataset who show similar ancestry to GBR reference individuals in R.
3. Using PLINK to identify a subset of individuals in the test dataset whose genotypes identify them as female and comparing to the reported sex values.
4. Evaluating the results of the sex-check and generating appropriate visualisations in R.
### Defining Ancestry Metrics
#### Soft Coding and Validation Methods
For the first part of the task, I developed a bash script in which the user could specify the paths for the reference data directory, test data directory, path to linkage disequilibrium (LD) ranges, output and results. I did this using optargs, shown below in *Code block 1* which allows the user to pass arguments to the script when invoking through the CLI. 

```bash
# Default values for ease of use for myself.
TEST=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results
REFERENCE=/home/c.c24082291/PLINK_Bioinfo/database/genotypes/all-1000g-phase3-chrall-mac5-v2/all-1000g-phase3-chrall-mac5-v2
HIGHLD=/home/c.c24082291/PLINK_Bioinfo/database/ranges/high-ld-b37.ranges
OUTDIR=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/tmp
RESULTS=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results
PLINK=""

# Parse command line arguments
while getopts ":t:r:h:o:s:p:" opt; do

    case ${opt} in

        t ) TEST=$OPTARG ;;

        r ) REFERENCE=$OPTARG ;;

        h ) HIGHLD=$OPTARG ;;

        o ) OUTDIR=$OPTARG ;;

        s ) RESULTS=$OPTARG ;;

		p ) PLINK=$OPTARG

        * ) usage ;;

    esac

done
```
*Code block 1*

To avoid confusion of how to pass arguments properly, I added a "usage" function, shown in *Code block 2* that would run if incorrect flags were passed on execution. This would output a help message which would inform a user how to structure their command. 
```bash
# Function to display help message

usage() {

    echo "Usage: $0 -t TEST -r REFERENCE -h HIGHLD -o OUTDIR -s RESULTS"

    echo "  -t TEST       Path to the test dataset"

    echo "  -r REFERENCE  Path to the reference dataset"

    echo "  -h HIGHLD     Path to the high LD ranges file"

    echo "  -o OUTDIR     Path to the output directory"

    echo "  -s RESULTS    Path to the results directory"

	echo "  -p PLINK      Path to PLINK"

    exit 1

}
```
*Code block 2*

Additionally I added simple data validation steps to ensure that the user provided arguments for each necessary flag shown in *Code block 3*.
```bash
# Check if all required arguments are provided

if [ -z "$TEST" ] || [ -z "$REFERENCE" ] || [ -z "$HIGHLD" ] || [ -z "$OUTDIR" ] || [ -z "$RESULTS" ]; then

    usage

fi
```
*Code block 3*
#### Main Script body

##### Extracting SNP IDs from both datasets

Now that the files necessary files and directories have been set, we want to start by identifying the SNP IDs in both the test and reference dataset. In the data used in this task that is stored in the 2nd column in the .bim file from both datasets, so using awk we can extract these and store separately. I also print the number of SNPs found in each to the console as an easy way to see if there are errors in reading the files e.g. 0 SNPs detected in either sample. The code for this is shown in *Code block 4*.

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

 In the next step, we use the `comm` command with the `-12` flag to compare the two files and keep only lines which are common to both test and reference data, shown in *Code block 5*. These common SNP IDs are stored in an `.extract` file. The number of common SNP IDs are printed to the console to inform the user of any unexpected absences of output. I could have added further validation to halt script execution if there were not any common SNPs.
 
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

The next step is to extract the data rows with the SNP information that are common to both reference and test datasets. This is different to the previous steps as we are extracting all of the related information not just the SNP ID. The common SNP data from each dataset is stored in a new `.bim` file for further processing as displayed in *Code block 6*.

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

The next step in the pipeline is to remove any SNPs with very low MAF. This is important as association signals from these SNPs are less informative as they occur in only a small number of individuals (Anderson, CA et al., 2010). The filtered dataset is then stored for further processing as demonstrated in *Code block 7.

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
These are considered ambiguous because without knowing the strand the same SNP could be reported as forward or reverse. For example a SNP reported as A/T on the forward strand could be T/A on the reverse strand. Without knowing the strand the true allele configuration is ambiguous and can lead to errors so it is important to exclude them from analysis. The code for this is shown below in *Code block 8*.

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

The next step is to limit the filtered SNPs to ancestry informative markers which are SNPs that have large frequency divergences among populations from different geographic regions and are can be used to estimate an individuals ancestry (Felkl, A.B et al., 2023). This is shown in *Code block 9* below. This also simplifies further analysis by focussing on a smaller set of markers which have the most importance.

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

Using plink we can create a set of SNPs that overlap with the ranges specified in our `ranges` file which are high LD regions as demonstrated in *Code block 10*. The SNPs are stored so they can be excluded in further analysis.

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

Next we use plink to generate a list of SNPs in a `.prune.in` file to prune from our datasets shown in *Code block 11*.

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

Using plink we can now extract the SNP information from the LD limited SNPs in the test and reference dataset in preparation to merge both shown in *Code block 12*.

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

The next step is to merge the two processed datasets. We begin by doing an initial merge attempt shown in *Code block 13*. If there are any strand alignment issues, the merge will fail and plink generates a list of SNPs which are problematic in a `.missnp` file.

```bash
# Merge datasets

plink \

--bfile "${OUTDIR}/reference_premerge" \

--bmerge "${OUTDIR}/test_premerge" \

--make-bed \

--out "${OUTDIR}/combined"
```
*Code block 13*

If the merge fails, the next step is to flip the strand of the SNPs with strand alignment issues and attempt the merge again as shown in *Code block 14*.

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

The penultimate step is to calculate principal components from the merged dataset using the command shown below in *Code block 15*. The output gives us an `.eigenvec` and `.eigenval` file with values from the matrix computation from the PCA.

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

The final step is to remove intermediatory files generated throughout the pipeline as shown below in *Code block 16*. This of course can be left out or modified to retain certain files if necessary. 

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

The first step after loading our libraries is to import our data into RStudio as a tibble using the `readr` library. As our data is headerless we need to define column names. Some columns need to have their datatype explicitly coerced such as FID and IID due to the front section of the samples having integer looking IDs which cause a type error when reading non-int-like rows included later in the data. R code for this is shown below in *Code block 17*.

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
#### Generating Scree Plot from `eigenvals`

After importing our data, we need to calculate explained variance for each component and plot them to create a Scree plot so that we can decide upon the number of principal components to retain. The code in *Code block 18* demonstrates how the dataset is created using the `eigenvalues`.

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

Next I used `ggplot2` and `plotly` libraries to plot the variance explained (plotted in blue) for each principal component as well as the cumulative variance (plotted in red) shown in *Code block 19*. I used `plotly` as it provides 'Show closest data on hover' which I find is helpful to explore values on the plot.

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

![Scree Plot](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/screeplot.png?raw=true)
*Fig 1, Scree Plot*

To determine the number of components to preserve we can for the "inflection point" or "elbow" in the Scree plot (Cangelosi, R. and Goriely, A., 2007). Looking at the output in *Fig 1*, that appears to be at the 3rd component so will only use data from the first 3 PCs.
#### Merge Population, Superpopulation and Eigenvectors

Next we want to visualise our principal component values to see if we can observe our unlabelled test data and how it is distributed compared to already know ancestral population/superpopulation groups.
The first step for this is to merge the population, superpopulation and eigenvector data. This will add ancestry population information to individuals in the reference dataset. Individuals in the test dataset will be given population and superpopulation values of "unknown" for the moment to distinguish them from the reference dataset. Code demonstrating this is shown below in *Code block 20*.

```r
pop_file <- read_table("E:/Bioinformatics/PCA_ASSESSMENT/all-1000g-phase3-chrall-mac5-v2.population", 
    col_names = c("FID", "IID", "Population"))

superpop_file <- read_table("E:/Bioinformatics/PCA_ASSESSMENT/all-1000g-phase3-chrall-mac5-v2.super-population", 
    col_names = c("FID", "IID", "Superpopulation"))
merged_pop_file <- full_join(pop_file, superpop_file)



merged_data <- eigenvec %>%
  left_join(merged_pop_file, by = "FID") %>%
  mutate(
    Population = if_else(is.na(Population), "unknown", Population),
    Superpopulation = if_else(is.na(Superpopulation), "unknown", Superpopulation),
  ) %>%
  select(-IID.y) %>%  # Remove the redundant IID.y column
  rename(IID = IID.x)  # Rename IID.x to IID
    

merged_data
```
*Code block 20*

#### Visualising PCs with 2-way Scatter Plots

Using `ggplot2` and our merged dataset we can now start plotting our components to see if we can observe trends as shown *Code block 21* below. I am using a colourblind friendly palette on the plots to facilitate visibility. I am also colouring by superpopulation for readability as using population would leave a very confusing plot with too many colours.

```r
ggplot(merged_data, aes(x = PC2, y = PC1, color = Superpopulation)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "PC2 vs PC1 (Colored by Superpopulation)",
    x = "PC2",
    y = "PC1"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
```
*Code block 21*

![PC1vsPC2](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/pc1vpc2.png?raw=true)
*Fig 2, PC2 vs PC1 plot*

As shown in *Fig 2* we can see that the majority of our test dataset individuals (coloured `unknown`) do not align to EAS or AFR ancestry based on the PC data. The code for these plots remains almost identical so I am not going to add the code for each to this report. The clustering in *Fig 3* and *Fig 4* also seem to indicate the test sample is made of majority EUR and AMR ancestry with a smaller number of SAS.

![pc1vs3](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/pc1vspc3.png?raw=true)
*Fig 3, PC1 vs PC3 plot*

![pc2vs3](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/pc2vspc3.png?raw=true)
*Fig 4, PC2 vs PC3 plot*
#### Calculating Mean and SD for each PC to Define Cluster Thresholds

As we want to identify GBR-like individuals in our test dataset, we want to extract the mean and standard deviation for each of our preserved principal components in order to define thresholds for defining if an individual in "GBR-like". This is shown below in *Code block 22*.

```r
gbr_like_pca_stats <- merged_data %>%
  filter(Population == "GBR") %>%
  select(starts_with("PC")) %>%
  summarise(across(everything(), list(mean = ~mean(.), sd = ~sd(.))))

gbr_like_pca_stats
```
*Code block 22*

#### Defining Thresholds

Finding the "correct" threshold to use is not an exact science and in the next steps I will demonstrate how I came to the number of standard deviations I used for defining my thresholds. As shown in *Code block 23*, I am using 2 standard deviations from the mean of the individuals defined as GBR population from the reference dataset.

```r
# Only use first 3 PCs as those explain most variance in data
# Within 2 sd
# k = num standard deviations from mean as threshold
k <- 2
thresholds <- tibble(
  PC1_mean_lower = gbr_like_pca_stats$PC1_mean - k * gbr_like_pca_stats$PC1_sd,
  PC1_mean_upper = gbr_like_pca_stats$PC1_mean + k * gbr_like_pca_stats$PC1_sd,
  PC2_mean_lower = gbr_like_pca_stats$PC2_mean - k * gbr_like_pca_stats$PC2_sd,
  PC2_mean_upper = gbr_like_pca_stats$PC2_mean + k * gbr_like_pca_stats$PC2_sd,
  PC3_mean_lower = gbr_like_pca_stats$PC3_mean - k * gbr_like_pca_stats$PC3_sd,
  PC3_mean_upper = gbr_like_pca_stats$PC3_mean + k * gbr_like_pca_stats$PC3_sd,
)
thresholds
```
*Code block 23*

#### Identifying Individuals who Fall Within Thresholds

The next step was to filter the dataset to find individuals whose PC values fall within the defined thresholds and store their cluster status in a new column in the `eigenvec` tibble as shown in *Code block 24* below.

```r
eigenvec <- eigenvec %>%
  rowwise() %>%
  mutate(
    cluster_status = ifelse(
      (PC1 >= thresholds$PC1_mean_lower & PC1 <= thresholds$PC1_mean_upper) &
      (PC2 >= thresholds$PC2_mean_lower & PC2 <= thresholds$PC2_mean_upper) &
      (PC3 >= thresholds$PC3_mean_lower & PC3 <= thresholds$PC3_mean_upper),
      "GBR-like", "Other"
    )
  )

# keep copy of all "gbr-like"
gbr_like_individuals <- eigenvec %>%
  filter(cluster_status == "GBR-like")

gbr_like_individuals
eigenvec
```
*Code block 24*

After filtering the data to see how many individuals across both datasets are included in the "GBR-like" cluster we have 1667 people. In the reference data there are 91 individuals labelled as being from the GBR population. Despite this being the case, when the clustering is ran with the 2 SD threshold, an additional 84 people from the reference dataset are identified as being "GBR-like". This mostly includes individuals labelled as CEU, which mean they are of northern European decent but the samples are from North America. This likely means that their ancestors emigrated from from Britain at some point in the recent past so this is not much of an issue. 
Increasing the threshold to 3 SDs from the mean increases the number of "mislabelled" individuals to 130 and begins to include an increasingly large number of samples from the IBS and TSI regions. This is why I decided to go with 2 SD from the mean as my threshold.
The code shown below in *Code block 25* demonstrates how the reference individuals were isolated and used to check for false positives.

```r
merged_data <- left_join(merged_data, select(eigenvec, FID, cluster_status), by = "FID")
ref_data <- merged_data %>%
  filter(Population != "unknown")
ref_data

pop_gbr_count <- ref_data %>%
  filter(Population == "GBR") %>%
  nrow()
pop_gbr_count
identified_gbr_count <- ref_data %>%
  filter(cluster_status == "GBR-like")
identified_gbr_count

```
*Code block 25*

Shown below in *fig 5* is the GBR-like cluster highlighted with the first 2 PCs. 

![gbrpc1vs2](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/gbr-highlighted.png?raw=true)
*Fig 5, GBR-like highlighted PC1vsPC2 plot*

```r
ggplot(eigenvec, aes(x = PC1, y = PC2, color = cluster_status)) +
  geom_point(alpha = 0.6, size = 3) +  # Scatter plot of individuals
  stat_ellipse(data = eigenvec %>% filter(cluster_status== "GBR-like"),
               aes(x=PC1, y=PC2),
               level = 0.95,
               color="blue", linewidth=1, linetype="dashed"
    ) +
  scale_color_manual(values = c("GBR-like" = "blue", "Other" = "gray")) +
  theme_minimal() +
  labs(title = "PCA Clustering with GBR 2-SD Threshold",
       x = "PC1", y = "PC2", color = "Cluster Status")
```
*Code block 26*

Code for *Fig 5* is shown in *Code block 26*. Similarly code for *Fig 6* is shown in *Code block 27*.

For once the data can be also be visualised in a 3D plot using `ggplot2` and `plotly` together to create an interactive plot which can be imbedded into the page shown as *Fig 6* (if working). Manipulation and zooming on this plot can better show the clustering in action.
<iframe src="https://rawcdn.githack.com/PlsNothingEmbarrassing/ancestry-pipeline/de50562300bae1833965ffd4dc2270b3855ed4fb/data/3dPCAplot.html" width="800" height="600"></iframe>
*Fig 6, 3D interactive PC plot highlighting GBR-individuals*

```r
library(rgl)
library(car)

fig <- plot_ly(eigenvec,
        x=~PC1, y=~PC2, z=~PC3,
        color = ~cluster_status,
        colors = c("GBR-like"= "blue", "Other" = "grey"),
        type = "scatter3d",
        mode="markers") %>%
  layout(title= "PCA plot",
         scene=list(
           xaxis=list(title="PC1"),
           yaxis=list(title="PC2"),
           zaxis=list(title="PC3")
         ))
htmlwidgets::saveWidget(fig, "3dPCAplot.html")
```
*Code block 27*

Finally we can isolate the test dataset individuals who are determined to be "GBR-like" and store their IDs in a `.keep` file as shown in *Code block 28*.

```r
merged_data_ref <- merged_data %>%
  filter(Population != "unknown")
merged_data_ref

test_gbr_like <- gbr_like_individuals %>%
  anti_join(merged_data_ref, by='FID')
test_gbr_like

test_gbr_keep <- test_gbr_like[1:2]
write.table(test_gbr_keep,
            file = "test_gbr.keep",
            col.names = F,
            row.names = F,
            quote = F,
            sep = " ")
```
*Code block 28*
### Identifying Females in subset

To identify females from our data subset we can use the PLINK `check-sex` method. Using bash we can create a script which we can pass our generated file.
The first section is some simple data validation methods to ensure the correct file is passed to the script as shown in *Code block 29*.

```bash
# Check if the user provided a keep file

if [ -z "$1" ]; then

    echo "Error: No .keep file specified."

    exit 1

fi 

# Assign the first argument to KEEP_FILE

KEEP_FILE=$1 

# Check if the keep file exists

if [ ! -f "$KEEP_FILE" ]; then

    echo "Error: The specified .keep file does not exist."

    exit 1

fi  

# Check if the keep file is empty

if [ ! -s "$KEEP_FILE" ]; then

    echo "Error: The specified .keep file is empty."

    exit 1

fi 

# Check if the keep file has the correct extension

if [[ "$KEEP_FILE" != *.keep ]]; then

    echo "Error: The specified file is not in the correct format."

    exit 1

fi
```
*Code block 29*

We then get the genome data associated with our subset individuals by filtering based on the `keep` file IDs as shown in *Code block 30*.
```bash
# Run plink to filter the dataset based on the keep file and create a new binary file
plink --bfile $TEST \

--keep "$KEEP_FILE" \

--make-bed \

--out "${OUTDIR}/test_gbr"
```
*Code block 30*

We then need to split the X chromosome pseudo-autosomal region off as they are sections which show higher heterozygosity in males than other regions which can make it harder to differentiate the genders. Code for this is shown in *Code block 31*.

```bash
# Split the X chromosome pseudo-autosomal region

plink --bfile "${OUTDIR}/test_gbr" \

--split-x b37 no-fail \

--make-bed \

--out "${OUTDIR}/test_gbr_split"
```
*Code block 31*

Next we can perform the sex check with the genome data then use awk to isolate only the individuals identified as female and store them in a separate file. Code shown in *Code block 32*.

```bash
# Run plink to check the sex of individuals in the filtered dataset

plink --bfile "${OUTDIR}/test_gbr_split" \

--check-sex \

--out "${OUTDIR}/test_gbr_sexcheck"  

# Extract individuals identified as female sex id =2 from the sex check results

awk '$4 == 2' "${OUTDIR}/test_gbr_sexcheck.sexcheck" > "${OUTDIR}/test_gbr_female.sexcheck"
```
*Code block 32*

The sex-check pipeline identified 694 genetic females and 723 males with 0 mismatches from reported sex.
### Evaluating Sex-Check Data

After gathering the sex-check data from our bash script, we can now import it into R Studio to evaluate and visualise elements of it.

The first step was to import and add column information as shown below in *Code block 33*.

```r
# Read in sex data for gbr test pop
sex_data <- read_table("test_gbr_sexcheck.sexcheck")
gbr_females <- read_table("test_gbr_female.sexcheck", col_names = c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F"))
females_keep <- gbr_females %>%
  select
write.table(females_keep,
            file = "test_females_gbr.keep",
            col.names = F,
            row.names = F,
            quote = F,
            sep = )
```
*Code block 33*

After the data was imported we can start making some visualisations. *Fig 7* displays the distributions across both sexes from the sex-check data.

![Sex Distrubtion](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/sex_distribution.png?raw=true)
*Fig 7, Sex Distribution*

Visualising X chromosome heterozygosity, stored as the F value in the plink sex-check data can be used to spot any clear outliers. There are not any clear outliers as shown below in *Fig 8*.

![X chromosome heterozygosity plot](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/X_heterozygosity.png?raw=true)
*Fig 8, X Chromosome Heterozygosity histogram*
## Data
### 1000 Genomes Project Data
The 1000 genomes project data is used as the reference dataset for my ancestry pipeline. The data is comprised of 2504 individuals from 26 populations representing a range of global populations.
#### Data Components
1. **Population files:** Population file contains individual ID (IID), family ID (FID), and Population label e.g. GBR, IBS.
2. **Superpopulation file:** Superpopulation file contains individual ID (IID), family ID (FID), and Superpopulation label e.g. EUR, SAS.
3. **Genotype Data:** Genotype data is stored in .bed, .bim and .fam files. Reported sex is stored in .fam file. Data contains SNP data along with other variant types.
4. **Ancestry Informative Markers (AIMs):** AIMs data contains list of SNPs that have been selected as Ancestry Informative Markers.
## Software
### PLINK

The PLINK version used in the bash scripts was PLINK 1.9 on the HAWK supercomputer.
### R Studio

The R code was written and run with R version 4.4.2 on R Studio.

The following libraries were used:
- `readr`
- `dplyr`
- `ggplot2`
- `plotly`
- `RcolorBrewer`
- `rgl`
- `car`

## Description of Workflow

Shown below in *Fig 9* are the steps taken for the ancestry analysis workflow.

![pca workflow diagram](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/pca_workflow.png?raw=true)
*Fig 9, Ancestry Analysis Workflow Diagram*

Shown below in *Fig 10* are the steps taken for the evaluating the ancestry data and finding individuals from the test dataset who fall into "GBR-like" ancestry.

![evaluate ancestry workflow diagram](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/evaluate_metrics.png?raw=true)
*Fig 10, Evaluating Ancestry Workflow Diagram*

*Fig 11* displayed below details the steps taken in to check the sex for the data subset.

![sex check workflow](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/sex_check_workflow.png?raw=true)
*Fig 11, Sex Check Workflow Diagram*

Shown in *Fig 12* is the workflow for the female evaluation.

![Evaluate female workflow](https://github.com/PlsNothingEmbarrassing/ancestry-pipeline/blob/main/data/evaluate_female.png?raw=true)
*Fig 12, Evaluate female Workflow Diagram*
## Discussion of Results

### Defining Ancestry Metrics Results

The result of the script defining ancestry metrics are the resultant eigen vector and eigen value files from the PCA which are then used in the second part of this assignment.

### Evaluating Results of PCA and Identifying GBR-like Individuals

There were 1492 GBR-like individuals identified from the test dataset which fell within our 2 SD cluster. Based on the PC visualisations this is not surprising as the majority of samples from the test dataset appeared to be relatively close to the European superpopulation values within our plotted PCs (shown in *Fig 2, 3 and 4*).

### Identifying Females in Subset and Evaluating

There were no abnormal results in the sex-check performed on the GBR-like individuals. The genomic sex of each individual matched the reported sex in the `.fam` file. Additionally after visualising the X chromosome heterozygosity, shown in *fig 8*, there are no individuals who appear wildly out of expected ranges.
## References

-  Anderson, C.A., Pettersson, F.H., Clarke, G.M., Cardon, L.R., Morris, A.P. and Zondervan, K.T., 2010. Data quality control in genetic case-control association studies. _Nature protocols_, _5_(9), pp.1564-1573.
- Byun, J., Han, Y., Gorlov, I.P., Busam, J.A., Seldin, M.F. and Amos, C.I., 2017. Ancestry inference using principal component analysis and spatial analysis: a distance-based analysis to account for population substructure. _BMC genomics_, _18_, pp.1-12.
- Cangelosi, R. and Goriely, A. (2007). Component retention in principal component analysis with application to cDNA microarray data. _Biology Direct_, 2(1), p.2. doi:https://doi.org/10.1186/1745-6150-2-2.
- Felkl, A.B., Avila, E., Gastaldo, A.Z., Lindholz, C.G., Dorn, M. and Alho, C.S. (2023). Ancestry resolution of South Brazilians by forensic 165 ancestry-informative SNPs panel. _Forensic Science International: Genetics_, [online] 64, p.102838. doi:https://doi.org/10.1016/j.fsigen.2023.102838.
-  Turner, S., Armstrong, L.L., Bradford, Y., Carlson, C.S., Crawford, D.C., Crenshaw, A.T., de Andrade, M., Doheny, K.F., Haines, J.L., Hayes, G. and Jarvik, G., 2011. Quality control procedures for genome‐wide association studies. _Current protocols in human genetics_, _68_(1), pp.1-19.
-  www.cog-genomics.org. (n.d.). _Basic statistics - PLINK 1.9_. [online] Available at: https://www.cog-genomics.org/plink/1.9/basic_stats.‌
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