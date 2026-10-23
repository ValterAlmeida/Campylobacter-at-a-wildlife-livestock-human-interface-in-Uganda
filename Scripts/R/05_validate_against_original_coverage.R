# =============================================================================
# Validation of the competitive-mapping abundances against the normalised
# coverage used in the original analysis (norm_jgi_cov_table.py output stored in
# bin_coverage_taxonomy.tsv). Pooled libraries only (UG_DNA_01 to UG_DNA_30).
#
# Output: Spearman correlation (overall and per pool) and a log-log scatter plot.
# Run after 03_MAG_relative_abundance.R.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2)
})

in_dir  <- "~/Downloads/MAG_relative_abundance"
out_dir <- file.path(in_dir, "R_output")
f_tax   <- "~/Downloads/R_coverage/bin_coverage_taxonomy.tsv"

cov_old <- read_tsv(f_tax, show_col_types = FALSE) %>%
  select(-classification) %>%
  pivot_longer(-bin_id, names_to = "col", values_to = "norm_cov_old") %>%
  mutate(sample = str_extract(col, "UG_DNA_\\d+")) %>%
  select(mag = bin_id, sample, norm_cov_old)

new <- read_tsv(file.path(out_dir, "pooled_mag_per_pool.tsv"), show_col_types = FALSE) %>%
  filter(filter == "all") %>% select(mag = unit, sample, norm_abund)

chk <- inner_join(cov_old, new, by = c("mag", "sample"))
message("MAG x pool pairs compared: ", nrow(chk))

rho <- cor(chk$norm_cov_old, chk$norm_abund, method = "spearman")
print(cor.test(chk$norm_cov_old, chk$norm_abund, method = "spearman", exact = FALSE)[c("estimate", "p.value")])
print(chk %>% group_by(sample) %>%
        summarise(rho = suppressWarnings(cor(norm_cov_old, norm_abund, method = "spearman"))) %>%
        summarise(median_rho = median(rho, na.rm = TRUE), min_rho = min(rho, na.rm = TRUE)))

# Zeros cannot sit on a log axis: place them at half the smallest positive value
floor_x <- min(chk$norm_cov_old[chk$norm_cov_old > 0]) / 2
floor_y <- min(chk$norm_abund[chk$norm_abund > 0]) / 2
chk <- chk %>% mutate(x = pmax(norm_cov_old, floor_x), y = pmax(norm_abund, floor_y))

p <- ggplot(chk, aes(x, y)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#2C5F8A") +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x = min(chk$x), y = max(chk$y), hjust = 0, vjust = 1,
           label = sprintf("Spearman ρ = %.2f\nn = %d MAG–pool pairs", rho, nrow(chk))) +
  labs(x = "Normalised coverage (original analysis)",
       y = "Normalised abundance (competitive mapping)") +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "validation_scatter_pooled.pdf"), p, width = 4, height = 4)
