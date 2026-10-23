# Clean, publication-ready Excel workbook of relative abundance (pooled + individual).
# Run after 03_MAG_relative_abundance.R. Reads only the two Supplementary_table_*.tsv files.
library(dplyr); library(readr); library(writexl)

out_dir <- "~/Downloads/MAG_relative_abundance/R_output"
xlsx    <- file.path(out_dir, "Supplementary_Table_relative_abundance.xlsx")

host_levels <- c("Human, healthcare", "Human, community", "Gorilla, habituated",
                 "Gorilla, non-habituated", "Cattle", "Goat")

tidy <- function(f, level) {
  read_tsv(file.path(out_dir, f), show_col_types = FALSE) %>%
    mutate(host_group = factor(host_group, host_levels)) %>%
    { stopifnot("unrecognised host_group label" = !anyNA(.$host_group)); . } %>%
    arrange(filter, species, host_group) %>%
    rename(`Host group` = host_group, Species = species,
           `Pools (n)` = pooled_n_pools,
           `Pooled mean (% of host-filtered reads)` = pooled_mean_pct,
           `Pooled min (%)` = pooled_min_pct, `Pooled max (%)` = pooled_max_pct,
           `Individual libraries (n)` = individual_n_libraries,
           `Prevalence (% of libraries with >=10 reads)` = individual_prevalence_pct,
           `Individual median (%)` = individual_median_pct,
           `Individual Q1 (%)` = individual_q25_pct, `Individual Q3 (%)` = individual_q75_pct,
           `Individual read-weighted (%)` = individual_aggregate_pct) %>%
    { if (level == "mag") rename(., MAG = mag, `MAG length (bp)` = length_bp) else . } %>%
    relocate(Species)
}
sheets <- list()
for (lv in c("species", "mag")) {
  f  <- if (lv == "mag") "Supplementary_table_relative_abundance_by_MAG.tsv"
        else "Supplementary_table_relative_abundance_by_species.tsv"
  d  <- tidy(f, lv)
  nm <- if (lv == "mag") "By MAG" else "By species"
  sheets[[paste0(nm, " (95% identity)")]] <- d %>% filter(filter == "id95") %>% select(-filter)
  sheets[[paste0(nm, " (all alignments)")]] <- d %>% filter(filter == "all") %>% select(-filter)
}
sheets[["Notes"]] <- data.frame(Notes = c(
  "Relative abundance of the 23 dereplicated Campylobacter MAGs, from competitive mapping of host-filtered reads.",
  "Pooled: 30 AGRF NovaSeq pools (5 per host group). Individual: 549 Orion MiSeq libraries with at least 10,000 read pairs (4 of 553 excluded).",
  "Abundance is the percentage of host-filtered reads assigned to the MAG or species (each read assigned to one MAG).",
  "Main sheets count alignments with at least 95% identity (NM/read length <= 0.05). 'All alignments' sheets are the sensitivity analysis.",
  "Species sums the reads of all MAGs assigned to one GTDB species; four gorilla MAGs without a species name are shown separately.",
  "Pooled values: mean, minimum and maximum across the pools of the host group.",
  "Individual values: median and interquartile range across libraries (mostly zero); read-weighted = total reads / total library size.",
  "Prevalence: libraries with at least 10 reads at the species or MAG level.",
  "Values are rounded to 3 significant figures."))

ord <- c("Notes", setdiff(names(sheets), "Notes"))
write_xlsx(sheets[ord], xlsx, format_headers = TRUE)
message("Written: ", normalizePath(xlsx))
