# =============================================================================
# Relative abundance of the 23 dereplicated Campylobacter MAGs
#   - pooled libraries     (UG_DNA_01 to UG_DNA_30, 5 pools per host group)
#   - individual libraries (553 ugNNN libraries, by host group)
#
# Input: the two long tables written by 02_map_reads_to_campy_mags.sh
#   mag_counts_long_individual.tsv   (sample, mag, mapped_reads, mapped_reads_id95, total_pairs)
#   mag_counts_long_pooled.tsv       (same columns)
#
# Definitions
#   library size (reads)  = 2 x total_pairs (every mate is counted when mapped)
#   % of reads            = 100 x mapped reads / library size
#   normalised abundance  = (mapped reads x 1e9) / (library size x MAG length bp)
#                           (the RPKM-like metric used in the manuscript methods)
#   filter "all"          = every mapped read
#   filter "id95"         = mapped reads with >= 95% identity to the MAG (NM/read length <= 0.05)
#   detected              = at least DETECT_MIN_READS reads (id95 filter, so low-level
#                           cross-recruitment from conserved regions is not counted)
# Each read is assigned to one MAG by Bowtie2, so summing reads over MAGs does not
# double count. Reads shared by near-identical MAGs are split arbitrarily between
# them, so the species-level tables are the safer ones to interpret.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# ---- 1. Paths and parameters (edit here) ------------------------------------
in_dir  <- "~/Downloads/MAG_relative_abundance"        # folder with the NeSI outputs
f_ind   <- file.path(in_dir, "mag_counts_long_individual.tsv")
f_pool  <- file.path(in_dir, "mag_counts_long_pooled.tsv")   # skipped if absent
f_len   <- file.path(in_dir, "campy_mags_23_lengths.tsv")
f_meta  <- "~/Downloads/Metagenomics_2026/Global_comparison/sample_metadata.tsv"
f_tax   <- "~/Downloads/R_coverage/bin_coverage_taxonomy.tsv"
out_dir <- file.path(in_dir, "R_output")

MIN_PAIRS_IND   <- 10000   # individual libraries smaller than this are excluded
DETECT_MIN_READS <- 10     # reads (id95) needed to call a MAG detected in a library

# Optional: readable names for species labels, e.g. c("bin.63.hh11" = "C. infans")
# applied to the species of the named MAG after the GTDB name is read
species_override <- c()

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Reference tables -----------------------------------------------------
len <- read_tsv(f_len, col_names = c("mag", "length_bp"), show_col_types = FALSE)

tax <- read_tsv(f_tax, show_col_types = FALSE) %>%
  select(mag = bin_id, classification) %>%
  mutate(species = str_remove(str_extract(classification, "s__[^;]*"), "^s__"),
         species = ifelse(is.na(species) | species == "", paste0("unnamed (", mag, ")"), species))
if (length(species_override) > 0) {
  tax$species[match(names(species_override), tax$mag)] <- unname(species_override)
}

message("MAG -> species used for the species-level tables (check this):")
print(as.data.frame(tax %>% left_join(len, by = "mag") %>% select(mag, species, length_bp)),
      row.names = FALSE)

pool_group <- function(n) {
  case_when(n <= 5  ~ "gorilla_wild",
            n <= 10 ~ "gorilla_habituated",
            n <= 15 ~ "human_healthcare",
            n <= 20 ~ "human_community",
            n <= 25 ~ "cattle",
            TRUE    ~ "goat")
}

# ---- 3. Read counts and build the long table ---------------------------------
read_counts <- function(f, dataset) {
  d <- read_tsv(f, show_col_types = FALSE) %>%
    left_join(len, by = "mag") %>%
    left_join(select(tax, mag, species), by = "mag") %>%
    mutate(lib_reads = 2 * total_pairs, dataset = dataset)
  stopifnot(!anyNA(d$length_bp), !anyNA(d$species))
  d
}

ind <- read_counts(f_ind, "individual")
meta <- read_tsv(f_meta, show_col_types = FALSE) %>% select(sample = sample_id, group)
ind <- left_join(ind, meta, by = "sample")

# Livestock libraries are not in sample_metadata.tsv. Assign them from the sample
# ID ranges: cattle ug401-ug484 (ug459 excluded; 83 libraries = pools 21-25) and
# goat ug501-ug573 (73 libraries = pools 26-30). Check these ranges against
# Supplementary Data 2 before you rely on them.
id_num <- as.integer(str_extract(ind$sample, "\\d+"))
ind$group <- case_when(!is.na(ind$group)            ~ ind$group,
                       id_num >= 401 & id_num <= 484 ~ "cattle",
                       id_num >= 501 & id_num <= 573 ~ "goat",
                       TRUE                          ~ NA_character_)
n_by_group <- ind %>% distinct(sample, group) %>% count(group)
message("Libraries per host group before the size filter (cattle should be 83, goat 73):")
print(n_by_group)

if (anyNA(ind$group)) {
  warning("Libraries without a host group: ",
          paste(unique(ind$sample[is.na(ind$group)]), collapse = ", "))
  ind <- filter(ind, !is.na(group))
}

n_all <- n_distinct(ind$sample)
small <- ind %>% distinct(sample, total_pairs) %>% filter(total_pairs < MIN_PAIRS_IND)
message(sprintf("Individual libraries: %d; excluded (< %d pairs): %d",
                n_all, MIN_PAIRS_IND, nrow(small)))
ind <- filter(ind, total_pairs >= MIN_PAIRS_IND)
message("Host groups (individual libraries kept):")
print(ind %>% distinct(sample, group) %>% count(group))

have_pool <- file.exists(f_pool)
if (have_pool) {
  pool <- read_counts(f_pool, "pooled") %>%
    mutate(group = pool_group(as.integer(str_extract(sample, "\\d+"))))
  message("Pools: ", n_distinct(pool$sample))
}

# ---- 4. Metrics at MAG or species level --------------------------------------
to_units <- function(d, level) {
  long <- bind_rows(
    d %>% mutate(filter = "all",  reads = mapped_reads),
    d %>% mutate(filter = "id95", reads = mapped_reads_id95))
  if (level == "mag") {
    out <- long %>% transmute(dataset, sample, group, filter, unit = mag,
                              lib_reads, total_pairs, reads, length_bp)
  } else {
    out <- long %>%
      group_by(dataset, sample, group, filter, unit = species, lib_reads, total_pairs) %>%
      summarise(reads = sum(reads), length_bp = mean(length_bp), .groups = "drop")
  }
  out %>% mutate(pct_reads = 100 * reads / lib_reads,
                 norm_abund = reads * 1e9 / (lib_reads * length_bp))
}

# individual: per host group (and all individuals together)
summarise_individual <- function(u) {
  bind_rows(u, mutate(u, group = "All individuals")) %>%
    group_by(filter, group, unit) %>%
    summarise(
      n_libraries   = n(),
      n_detected    = sum(reads >= DETECT_MIN_READS),
      prevalence_pct = 100 * n_detected / n_libraries,
      median_pct    = median(pct_reads),
      q25_pct       = quantile(pct_reads, 0.25),
      q75_pct       = quantile(pct_reads, 0.75),
      max_pct       = max(pct_reads),
      aggregate_pct = 100 * sum(reads) / sum(lib_reads),
      median_pct_detected = ifelse(n_detected > 0,
                                   median(pct_reads[reads >= DETECT_MIN_READS]), NA_real_),
      median_norm_abund = median(norm_abund),
      .groups = "drop")
}

# pooled: per pool, then per host group (5 pools) and over all pools
summarise_pooled <- function(u) {
  by_group <- bind_rows(u, mutate(u, group = "All pools")) %>%
    group_by(filter, group, unit) %>%
    summarise(
      n_pools       = n(),
      mean_pct      = mean(pct_reads),
      min_pct       = min(pct_reads),
      max_pct       = max(pct_reads),
      aggregate_pct = 100 * sum(reads) / sum(lib_reads),
      mean_norm_abund = mean(norm_abund),
      .groups = "drop")
  by_group
}

for (lev in c("mag", "species")) {
  u_ind <- to_units(ind, lev)
  s_ind <- summarise_individual(u_ind)
  write_tsv(s_ind, file.path(out_dir, paste0("individual_", lev, "_by_group.tsv")))
  write_tsv(u_ind,  file.path(out_dir, paste0("individual_", lev, "_per_library.tsv")))

  if (have_pool) {
    u_pool <- to_units(pool, lev)
    write_tsv(u_pool, file.path(out_dir, paste0("pooled_", lev, "_per_pool.tsv")))
    write_tsv(summarise_pooled(u_pool), file.path(out_dir, paste0("pooled_", lev, "_by_group.tsv")))
  }
}

# ---- 5. Total Campylobacter share (all 23 MAGs) ------------------------------
tot <- function(d) {
  bind_rows(d %>% mutate(filter = "all",  reads = mapped_reads),
            d %>% mutate(filter = "id95", reads = mapped_reads_id95)) %>%
    group_by(dataset, sample, group, filter, lib_reads) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    mutate(pct_reads = 100 * reads / lib_reads)
}
tot_all <- bind_rows(tot(ind), if (have_pool) tot(pool))
write_tsv(tot_all, file.path(out_dir, "total_campylobacter_pct_per_library.tsv"))
message("Total % of host-filtered reads mapping to the 23 MAGs (id95), by dataset and group:")
print(tot_all %>% filter(filter == "id95") %>%
        group_by(dataset, group) %>%
        summarise(n = n(), median = median(pct_reads), q25 = quantile(pct_reads, .25),
                  q75 = quantile(pct_reads, .75), max = max(pct_reads), .groups = "drop"),
      n = 40)

# ---- 6. Diagnostic: how much recruitment is below 95% identity ---------------
diag <- bind_rows(ind, if (have_pool) pool) %>%
  group_by(dataset, mag, species) %>%
  summarise(reads_all = sum(mapped_reads), reads_id95 = sum(mapped_reads_id95),
            frac_id95 = reads_id95 / reads_all, .groups = "drop") %>%
  arrange(dataset, frac_id95)
write_tsv(diag, file.path(out_dir, "diagnostic_fraction_reads_id95.tsv"))
message("MAGs with the lowest fraction of reads at >= 95% identity (individual libraries):")
print(diag %>% filter(dataset == "individual") %>% head(8))

# ---- 7. Host labels (shared by the figures and the supplementary tables) ----
# Pooled "gorilla_wild" is treated as the non-habituated gorilla group and pooled
# "human_healthcare/community" as "ugandan_healthcare/community". Check this
# against Supplementary Data 1 and 2.
host_ind  <- c(gorilla_nonhabituated = "Gorilla, non-habituated", gorilla_habituated = "Gorilla, habituated",
               ugandan_healthcare = "Human, healthcare", ugandan_community = "Human, community",
               cattle = "Cattle", goat = "Goat")
host_pool <- c(gorilla_wild = "Gorilla, non-habituated", gorilla_habituated = "Gorilla, habituated",
               human_healthcare = "Human, healthcare", human_community = "Human, community",
               cattle = "Cattle", goat = "Goat")
host_levels <- c("Human, healthcare", "Human, community", "Gorilla, habituated",
                 "Gorilla, non-habituated", "Cattle", "Goat")
sp_first <- c("Campylobacter infans", "Campylobacter_D sp900539255",
              "Campylobacter ovis", "Campylobacter sp017646085")
sp_levels <- c(intersect(sp_first, tax$species),
               setdiff(sort(unique(tax$species)), sp_first))

# ---- 8. Figures (supplementary) ----------------------------------------------
heat_theme <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(4, "mm"))

u_sp_ind <- to_units(ind, "species")
s_sp_ind <- summarise_individual(u_sp_ind) %>%
  filter(filter == "id95", group != "All individuals") %>%
  mutate(host = factor(unname(host_ind[group]), levels = host_levels),
         unit = factor(unit, levels = rev(sp_levels)))

p_ind <- ggplot(s_sp_ind, aes(host, unit, fill = prevalence_pct)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(prevalence_pct), colour = prevalence_pct > 50), size = 2.2,
            show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_fill_viridis_c(name = "Prevalence (%)", limits = c(0, 100)) +
  labs(x = NULL, y = NULL,
       title = "Individual libraries: % of libraries with the lineage detected") +
  heat_theme
ggsave(file.path(out_dir, "heatmap_individual_species.pdf"), p_ind, width = 150, height = 90,
       units = "mm", device = "pdf")

if (have_pool) {
  u_sp_pool <- to_units(pool, "species")
  s_sp_pool <- summarise_pooled(u_sp_pool) %>%
    filter(filter == "id95", group != "All pools") %>%
    mutate(host = factor(unname(host_pool[group]), levels = host_levels),
           unit = factor(unit, levels = rev(sp_levels)),
           mean_pct_plot = ifelse(mean_pct > 0, mean_pct, NA_real_),
           lab = ifelse(mean_pct < 0.001, "<0.001", sprintf("%.3f", mean_pct)),
           light = log10(mean_pct_plot) > mean(range(log10(mean_pct_plot), na.rm = TRUE)))
  p_pool <- ggplot(s_sp_pool, aes(host, unit, fill = mean_pct_plot)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = lab, colour = light), size = 2.2, na.rm = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "white"), na.value = "black") +
    scale_fill_viridis_c(name = "Mean % of reads", trans = "log10", na.value = "grey90") +
    labs(x = NULL, y = NULL, title = "Pooled libraries: mean % of host-filtered reads") +
    heat_theme
  ggsave(file.path(out_dir, "heatmap_pooled_species.pdf"), p_pool, width = 150, height = 90,
         units = "mm", device = "pdf")
}

# ---- 9. Supplementary table: pooled and individual side by side ---------------
make_supp_table <- function(level) {
  si <- summarise_individual(to_units(ind, level)) %>%
    filter(group != "All individuals") %>%
    transmute(filter, host_group = unname(host_ind[group]), unit,
              individual_n_libraries = n_libraries,
              individual_prevalence_pct = prevalence_pct,
              individual_median_pct = median_pct,
              individual_q25_pct = q25_pct, individual_q75_pct = q75_pct,
              individual_aggregate_pct = aggregate_pct)
  if (!have_pool) return(si)
  sp <- summarise_pooled(to_units(pool, level)) %>%
    filter(group != "All pools") %>%
    transmute(filter, host_group = unname(host_pool[group]), unit,
              pooled_n_pools = n_pools, pooled_mean_pct = mean_pct,
              pooled_min_pct = min_pct, pooled_max_pct = max_pct)
  full_join(sp, si, by = c("filter", "host_group", "unit"))
}

supp_mag <- make_supp_table("mag") %>%
  left_join(select(tax, mag, species), by = c("unit" = "mag")) %>%
  left_join(len, by = c("unit" = "mag")) %>%
  rename(mag = unit) %>%
  relocate(species, mag, length_bp) %>%
  arrange(filter, species, mag, host_group) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
supp_sp <- make_supp_table("species") %>%
  rename(species = unit) %>%
  arrange(filter, species, host_group) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
write_tsv(supp_mag, file.path(out_dir, "Supplementary_table_relative_abundance_by_MAG.tsv"))
write_tsv(supp_sp,  file.path(out_dir, "Supplementary_table_relative_abundance_by_species.tsv"))

message("Done. Tables and figures written to: ", normalizePath(out_dir))
