# Run from project directory
library(cancereffectsizeR)
library(data.table)

# Version 2.8 used (and Version 3.0 many include changes that will break this script)
stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('2.8.0') && 
          packageVersion('cancereffectsizeR') < as.package_version('3.0.0'))

# Genes reported in Brady figures 2A (B-ALL) and 2B (T-ALL)
genes_of_interest_B = read.table('source_data/genes_of_interest_B.txt', header = F, comment.char = '#')[[1]]
genes_of_interest_T = read.table('source_data/genes_of_interest_T.txt', header = F, comment.char = '#')[[1]]

# Read in prepped MAF data
maf = preload_maf(maf = "brady_filtered_hg38.maf", refset = "ces.refset.hg38", keep_extra_columns = 'Lineage')

# Verify that lineage (B vs. T) looks okay. A few rows have missing values, but they can be ignored.
lineage_info = unique(maf[! is.na(Lineage) & Lineage != '', .(Unique_Patient_Identifier, Lineage)])
stopifnot(lineage_info[, .N] == uniqueN(lineage_info$Unique_Patient_Identifier),
          uniqueN(maf$Unique_Patient_Identifier) == lineage_info[, .N],
          all(lineage_info$Lineage %in% c('T', 'B')))

# Create cancereffectsizeR analysis and load data
cesa = CESAnalysis(refset = "ces.refset.hg38")
cesa =  load_maf(cesa = cesa, maf = maf)
cesa = load_sample_data(cesa, lineage_info)

t_samples = cesa$samples[Lineage == "T", Unique_Patient_Identifier]
b_samples = cesa$samples[Lineage == "B", Unique_Patient_Identifier]

# Infer trinucleotide-context-specific relative rates of SNV mutation from a mutational signature analysis
# Leaving out signatures suggested absent in ALL by COSMIC.
# Running separately for T and B cell types.
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "ALL", treatment_naive = TRUE)
cesa <- trinuc_mutation_rates(
  cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.2,
  signature_exclusions = signature_exclusions,
  samples = t_samples
)
cesa <- trinuc_mutation_rates(
  cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.2,
  signature_exclusions = signature_exclusions,
  samples = b_samples
)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$lymphoma,
                            samples = t_samples)

cesa = gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$lymphoma,
                            samples = b_samples)

# Get variant counts by lineage and merge in gene annotation
all_counts = variant_counts(cesa, by = 'Lineage')
all_counts[cesa$variants, gene := gene, on = 'variant_id']

# Take all genes with at least 5 mutations
# And yes, we're going to leave out genes of interest with fewer mutations.
top_B_genes = all_counts[, sum(B_prevalence), by = 'gene'][V1 > 4, gene]
variants_to_use_B = all_counts[B_prevalence > 0][gene %in% top_B_genes, variant_id]
b_comp = define_compound_variants(cesa = cesa, variant_table = select_variants(cesa, variant_ids = variants_to_use_B),
                         by = "gene", merge_distance = Inf)
cesa = ces_variant(cesa, variants = b_comp, samples = b_samples, run_name = "comp_b")


top_T_genes = all_counts[, sum(T_prevalence), by = 'gene'][V1 > 4, gene]
variants_to_use_T = all_counts[T_prevalence > 0][gene %in% top_T_genes, variant_id]
t_comp = define_compound_variants(cesa = cesa, 
                                  variant_table = select_variants(cesa, variant_ids = variants_to_use_T),
                                  by = "gene", merge_distance = Inf)
cesa = ces_variant(cesa, variants = t_comp, samples = t_samples, run_name = "comp_t")
save_cesa(cesa, 'cesa.rds')

