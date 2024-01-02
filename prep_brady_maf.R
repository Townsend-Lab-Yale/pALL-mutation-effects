library(cancereffectsizeR)
library(ces.refset.hg38)
library(data.table)

# Get file for liftOver
chain_file_gz = 'hg19ToHg38.over.chain.gz' 
chain_file = 'hg19ToHg38.over.chain'
chain_url = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'

if(! file.exists(chain_file)) {
  download.file(chain_url, chain_file_gz)
  writeLines(readLines(chain_file_gz), chain_file)
}

# Read in Table S3 (sheet 3 of Excel file), skipping first lines until column headers
dt = as.data.table(readxl::read_excel('source_data/41588_2022_1159_MOESM4_ESM.xlsx', sheet = 3, skip = 3))

# First row after header is column definitions rather than data, so remove it
dt = dt[2:.N]

# Change column names to match MAF format
setnames(dt, old = c('Chr', 'Pos', 'Chr_Allele', 'Alternative_Allele', 'PatientID'),
         new = c('Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele', 'Unique_Patient_Identifier'))

maf = preload_maf(maf = dt, refset = 'ces.refset.hg38', chain_file = chain_file,
                  coverage_intervals_to_check = ces.refset.hg38$default_exome, keep_extra_columns = TRUE)

# Few problems, nothing concerning.
## problem                          num_records
## reference_mismatch               3
## duplicate_record_after_liftOver  2
## invalid_record                   1
## failed_liftOver                  31

# Remove problematic recrods.
maf = maf[is.na(problem)]

# Check how many SNV records are in or near coverage
snvs = maf[variant_type == 'snv']
coverage_stats =  snvs[, .(covered = mean(dist_to_coverage_intervals == 0), 
                      within_100 = mean(dist_to_coverage_intervals <= 100),
                      within_1000 = mean(dist_to_coverage_intervals <= 1000),
                      within_10k = mean(dist_to_coverage_intervals <= 10000))]
coverage_stats
# covered within_100 within_1000 within_10k
# 1: 0.8973903  0.9506393   0.9556918  0.9654233


# Most records are covered and another 5% are nearly covered (it's common for sequencing to pick up
# variants within 100bp of targets). Then, we pick up pick up records when expanding from 100bp to
# 1000bp out-of-coverage, suggesting remaining records are miscalled.
# Technically, it's possible that our exome definitions match theirs extremely well for the genes in our
# definitions (in that the 0/100/1000 pattern is just what you'd expect) and that they also
# incorporate additional coding regions that we don't use, but this seems unlikely.

# One more check: See what fraction of calls are in repetitive regions of the genome by distance to coverage.
rep_frac_by_dist = c(
  '<100' = snvs[dist_to_coverage_intervals <= 100, mean(repetitive_region)],
  '101-1000' = snvs[dist_to_coverage_intervals %in% 101:1000, mean(repetitive_region)],
  '1001-10000' = snvs[dist_to_coverage_intervals %in% 1001:10000, mean(repetitive_region)],
  '>10000' = snvs[dist_to_coverage_intervals > 10000, mean(repetitive_region)]
)
rep_frac_by_dist
# <100   101-1000 1001-10000     >10000 
# 0.02685153 0.10000000 0.18284424 0.16010165 

# Fraction of calls in repetitive regions goes up when you're >100 out-of-coverage,
# which is further evidence that the far out-of-coverage calls are low quality.

# Therefore, we'll remove records >100 bp out-of-coverage.
maf = maf[dist_to_coverage_intervals <= 100]

# We'll also do our standard filtering of removing variants that could be germline (rather than somatic)
# and removing all repetitive regions calls except for those at COSMIC-annotated sites. (COSMIC annotations
# indicates known or suspected cancer relevance, and annotated variants are assigned into tiers 1-3.)
maf = maf[germline_variant_site == FALSE]
maf = maf[repetitive_region == FALSE | cosmic_site_tier %in% 1:3]

# We kept ~91% of calls.
(filtering_summary = (c(kept = maf[, .N], original = dt[, .N])))
# kept original 
# 45037    49423 


fwrite(maf, 'brady_filtered_hg38.maf', sep = "\t")


