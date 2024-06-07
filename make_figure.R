library(cancereffectsizeR)
library(data.table)
library(ggplot2)

# Load the saved CESAnalysis that was generated in pALL_analysis.R 
cesa = load_cesa('cesa.rds')

# Genes reported in Brady figures 2A (B-ALL) and 2B (T-ALL)
genes_of_interest_B = fread('source_data/genes_of_interest_B.txt', header = F)[[1]]
genes_of_interest_T = fread('source_data/genes_of_interest_T.txt', header = F)[[1]]

b_results = cesa$selection$comp_b
b_results[, gene := gsub('\\.1', '', variant_name)] # get gene name from compound variant ID

# A few variants may appear in <4 samples because the variant counting above was not by sample; that is,
# some samples could have more than one of the variants, bringing the sample count below our filter.
# So, we just remove any variants with too few samples here.
b_results = b_results[included_with_variant > 4] # yes, dropping genes of interest with few mutations

num_B_new = 50 - b_results[gene %in% genes_of_interest_B, .N]
for_b_plot = rbind(b_results[gene %in% genes_of_interest_B],
                   b_results[order(-selection_intensity)][! gene %in% genes_of_interest_B][1:num_B_new])

for_b_plot = for_b_plot[order(-selection_intensity)] # sort ourselves to ensure that styling matches

# Bold the previously unreported genes. Note that ggplot doesn't like this way of styling labels
# and gives warnings, but as of ggplot v3.4.3, it does work.
for_b_plot[, to_style := 'bold']
for_b_plot[gene %in% genes_of_interest_B, to_style := 'plain']

# Repeat same process with T-cell group
t_results = cesa$selection$comp_t
t_results[, gene := gsub('\\.1', '', variant_name)] # get gene name from compound variant ID
t_results = t_results[included_with_variant > 4]

num_T_new = 50 - t_results[gene %in% genes_of_interest_T, .N]
for_t_plot = rbind(t_results[gene %in% genes_of_interest_T],
                   t_results[order(-selection_intensity)][! gene %in% genes_of_interest_T][1:num_T_new])

for_t_plot = for_t_plot[order(-selection_intensity)] # sort ourselves to ensure that styling matches
for_t_plot[, to_style := 'bold']
for_t_plot[gene %in% genes_of_interest_T, to_style := 'plain']

b_plot = plot_effects(for_b_plot, prevalence_method = 'both', y_label = 'gene', y_title = 'Gene',
                      topn = NULL, order_by_effect = FALSE, x_title = 'B-ALL cancer effect',
                      legend_size_name = 'Samples with mutated gene\n(Percent of B-ALL cohort)',
                      color_by = 'darkseagreen3', legend.position = c(.75, .20)) + 
  theme(axis.text.y = element_text(face = rev(for_b_plot$to_style)))

t_plot = plot_effects(for_t_plot, prevalence_method = 'both', y_label = 'gene', y_title = 'Gene',
                      topn = NULL, order_by_effect = FALSE, x_title = 'T-ALL cancer effect',
                      color_by = 'khaki', legend_size_name = 'Samples with mutated gene\n(Percent of T-ALL cohort)',
                      legend.position = c(.75, .20)) +
  theme(axis.text.y = element_text(face = rev(for_t_plot$to_style)))


## Previous submission used a combined plot:
# combined = cowplot::plot_grid(t_plot, b_plot, labels = 'AUTO', label_x = .03, label_y = .915, hjust = 0, vjust = 0)
# ggsave(plot = combined, filename = 'ALL_effects_plot.png', width = 1325*2.5, height = 750*2.5, 
#        units = 'px', dpi = 'retina')

# Save separate plots
ggsave(plot = b_plot, filename = 'ALL_effects_B-cell.png', width = 550 * 2.5, height = 800*2.5, units = 'px', dpi = 'retina')
ggsave(plot = t_plot, filename = 'ALL_effects_T-cell.png', width = 550 * 2.5, height = 800*2.5, units = 'px', dpi = 'retina')

