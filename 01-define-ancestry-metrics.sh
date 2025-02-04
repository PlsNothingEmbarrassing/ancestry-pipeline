#!/bin/bash
#SBATCH --account=scw1557
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=define-ancestry-metrics

#cd $SLURM_SUBMIT_DIR

# Default values for variables
TEST=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results
REFERENCE=/home/c.c24082291/PLINK_Bioinfo/database/genotypes/all-1000g-phase3-chrall-mac5-v2/all-1000g-phase3-chrall-mac5-v2
HIGHLD=/home/c.c24082291/PLINK_Bioinfo/database/ranges/high-ld-b37.ranges
OUTDIR=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/tmp
RESULTS=/home/c.c24082291/PLINK_Bioinfo/project/PCA_ASSESSMENT/results

# Function to display help message
usage() {
    echo "Usage: $0 -t TEST -r REFERENCE -h HIGHLD -o OUTDIR -s RESULTS"
    echo "  -t TEST       Path to the test dataset"
    echo "  -r REFERENCE  Path to the reference dataset"
    echo "  -h HIGHLD     Path to the high LD ranges file"
    echo "  -o OUTDIR     Path to the output directory"
    echo "  -s RESULTS    Path to the results directory"
    exit 1
}

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

# Check if all required arguments are provided
if [ -z "$TEST" ] || [ -z "$REFERENCE" ] || [ -z "$HIGHLD" ] || [ -z "$OUTDIR" ] || [ -z "$RESULTS" ]; then
    usage
fi


module load plink/1.9
# limit to snps in both datasets
echo "awking snps"
awk '{print $2}' "${REFERENCE}.bim" | sort > "${OUTDIR}/reference.snps"
awk '{print $2}' "${TEST}.bim" | sort > "${OUTDIR}/test.snps"

echo "Number of snps in reference dataset: $(wc -l < "${OUTDIR}/reference.snps")"
echo "Number of snps in test dataset: $(wc -l < "${OUTDIR}/test.snps")"

# Use comm to identify common snps
# Use -12 flag to suppress lines unique to each file
comm -12  \
"${OUTDIR}/reference.snps" \
"${OUTDIR}/test.snps" \
> "${OUTDIR}/common_to_both.extract"

echo "Number of snps in common: $(wc -l < "${OUTDIR}/common_to_both.extract")"

# limit to commons snps
# Use plink to extract common snps from ref dataset
# Use --extract flag to specify snps to keep
plink \
--bfile "${REFERENCE}" \
--extract "${OUTDIR}/common_to_both.extract" \
--make-bed \
--out "${OUTDIR}/reference_common"



# Use plink to extract common snps from in the new ref dataset
# Use maf flag to specify minor allele frequency
plink \
--bfile "${OUTDIR}/reference_common" \
--maf 0.05 \
--make-bed \
--out "${OUTDIR}/reference_common_maf5"
echo "Number of snps in common with MAF > 0.05: $(wc -l < "${OUTDIR}/reference_common_maf5.bim")"

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

# Limit to ancestry informative markers
# Use plink2 to extract ancestry informative markers from the ref dataset
plink \
--bfile "${OUTDIR}/reference_common_maf5_noWS" \
--extract "${REFERENCE}.aims" \
--make-bed \
--out "${OUTDIR}/reference_common_maf5_noWS_aims"

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

# LD prune
# Use plink to generate a LD pruned list of snps .prune.in
plink \
--bfile "${OUTDIR}/reference_common_maf5_noWS_aims_noLD" \
--indep-pairwise 50 10 0.1 \
--out "${OUTDIR}/ld_independent"
# Specify the parameters for LD pruning (window size, step size, r^2 threshold)

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

# Merge datasets
plink \
--bfile "${OUTDIR}/reference_premerge" \
--bmerge "${OUTDIR}/test_premerge" \
--make-bed \
--out "${OUTDIR}/combined"
# Use plink to update the reference dataset with the merged dataset
plink \
--bfile "${OUTDIR}/reference_premerge" \
--flip "${OUTDIR}/combined-merge.missnp" \
--make-bed \
--out "${OUTDIR}/reference_flipped"

# Use plink to merge the two ld limited datasets including flipped
plink \
--bfile "${OUTDIR}/reference_flipped" \
--bmerge "${OUTDIR}/test_premerge" \
--make-bed \
--out "${OUTDIR}/combined"

# Calculate principal components
# Use plink to calculate the principal components
plink \
--bfile "${OUTDIR}/combined" \
--pca \
--out "${RESULTS}/combined"

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
