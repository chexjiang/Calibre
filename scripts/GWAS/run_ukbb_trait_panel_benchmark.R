#!/usr/bin/env Rscript

# UKBB SCENT-style sensitivity analysis.
#
# This is a new presentation layer. It reuses the completed UKBB per-trait
# overlap tables and leaves the previous all-94 results untouched.
#
# Outputs:
#   primary: biology-matched panels, mean log2(1 + enrichment), trait bootstrap
#   secondary: all-94 traits, balanced across predefined trait domains
#   supplementary: all-94 traits, ordinary trait-level summary

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
out <- file.path(gwas, "SCENT_style", "UKBB")
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
for (d in c("config", "tables", "figures")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

n_boot <- if (length(args) >= 2L) as.integer(args[[2]]) else 1000L
seed <- if (length(args) >= 3L) as.integer(args[[3]]) else 20261002L

method_colors <- c(
  "EpiLink-Poisson" = "#D65F5F",
  "EpiLink-NB" = "#E39C37",
  "scMultiMap" = "#009E73",
  "Signac (corrected)" = "#7F7F7F",
  "Signac (uncorrected)" = "#2b8cbe"
)
methods <- names(method_colors)

canonical_method <- function(x) {
  y <- trimws(as.character(x))
  y[grepl("^EpiLink-NB", y, ignore.case = TRUE)] <- "EpiLink-NB"
  y[grepl("^EpiLink-Poisson", y, ignore.case = TRUE)] <- "EpiLink-Poisson"
  y[grepl("^scMultiMap$", y, ignore.case = TRUE)] <- "scMultiMap"
  unc <- grepl("^Signac.*uncorrected", y, ignore.case = TRUE)
  y[unc] <- "Signac (uncorrected)"
  y[!unc & grepl("^Signac.*corrected", y, ignore.case = TRUE)] <- "Signac (corrected)"
  y
}

read_per_trait <- function(path) {
  if (!file.exists(path)) stop("Missing per-trait UKBB result: ", path)
  x <- fread(path)
  x[, method_label := canonical_method(method_label)]
  x <- x[method_label %in% methods]
  x[, celltype := as.character(celltype_scope)]
  x[celltype == "ALL_REQUESTED", celltype := "ALL"]
  x[, view := ifelse(celltype == "ALL", "pooled", "celltype")]
  x[, method_label := factor(method_label, levels = methods)]
  x[, source_file := path]
  x
}

pooled <- read_per_trait(file.path(source_results$pooled_results,
                                  "scent_sakaue_per_trait.tsv.gz"))
celltype <- read_per_trait(file.path(source_results$celltype_results,
                                    "scent_sakaue_per_trait.tsv.gz"))

# Current benchmark scope: retain the original Brain dataset, add Brain2_ALL,
# Brain2_AD and Pituitary, and exclude NEAT/Brain2_CTRL from the new report.
pooled <- pooled[dataset != "NEAT"]
celltype <- celltype[dataset != "NEAT"]
pooled <- pooled[dataset != "Brain2_CTRL"]
celltype <- celltype[dataset != "Brain2_CTRL"]

uk <- rbindlist(list(pooled, celltype), fill = TRUE)
uk[, score := fifelse(is.finite(enrichment) & enrichment >= 0,
                      log2(1 + enrichment), NA_real_)]
uk[, dataset_celltype := ifelse(view == "pooled", dataset,
                                paste(dataset, celltype, sep = " | "))]

# Fixed, pre-specified biology-matched panels. These are written to disk below
# so the trait choice is auditable and cannot silently change with the result.
panels <- list(
  immune = c("AID_Combined", "Asthma", "IBD", "T2D", "CRP", "WBC", "Lym",
             "Neutro", "Mono", "Baso", "Eosino"),
  brain_AD = c("Alzheimer_LTFH"),
  endocrine = c("IGF1", "Testosterone", "SHBG", "Age_at_Menarche", "Age_at_Menopause"),
  stromal_exploratory = c("AID_Combined", "IBD", "CRP")
)

panel_map <- function(dataset, celltype) {
  if (dataset %in% c("Brain", "Brain2_ALL", "Brain2_AD")) return("brain_AD")
  if (dataset == "Pituitary") return("endocrine")
  if (dataset == "Arthritis" && celltype %in% c("Endothelial", "Fibroblast"))
    return("stromal_exploratory")
  "immune"
}

uk[, matched_panel := mapply(panel_map, dataset, celltype)]
panel_rows <- rbindlist(lapply(names(panels), function(p) {
  data.table(panel = p, trait = panels[[p]], panel_role = "primary")
}), fill = TRUE)
fwrite(panel_rows, file.path(out, "config", "biology_matched_trait_panels.tsv"), sep = "\t")

# Trait domains for the secondary all-94 domain-balanced score. The mapping is
# intentionally explicit; unlisted traits remain in "other" and are still
# included, so no trait is silently dropped from the all-94 sensitivity result.
domain_map <- list(
  immune_inflammatory = c("AID_Combined", "Asthma", "CRP", "IBD", "T2D", "T2D_BMI"),
  hematologic = c("Baso", "Eosino", "Hb", "Ht", "Lym", "MCH", "MCHC", "MCV", "Mono", "Neutro", "Plt", "RBC", "WBC"),
  cardiovascular = c("AFib", "CAD", "DBP", "DVT", "MAP", "PP", "SBP"),
  metabolic = c("AG", "Alb", "ALP", "ALT", "AST", "BFP", "BMI", "BW", "Ca", "eGFR", "eGFRcys", "Glucose", "HbA1c", "HDLC", "LDLC", "LipoA", "TBil", "TC", "TG", "TP", "UA", "Urea", "VitD", "WHRadjBMI"),
  endocrine_reproductive = c("Age_at_Menarche", "Age_at_Menopause", "IGF1", "SHBG", "Testosterone", "Hypothyroidism"),
  neuropsychiatric = c("Alzheimer_LTFH", "Depression_GP", "FedUp_Feelings", "Guilty_Feelings", "Insomnia", "Irritability", "Loneliness", "Migraine_Self", "Miserableness", "Mood_Swings", "Morning_Person", "Nervous_Feelings", "Neuroticism", "Risk_Taking", "Sensitivity", "Suffer_from_Nerves", "Tense", "Worrier", "Worry_Too_Long"),
  cancer = c("BrC", "CRC", "PrC"),
  behavior = c("Smoking_CPD", "Smoking_Ever_Never"),
  other = c("Balding_Type4", "Blood_Clot_Lung", "Cholelithiasis", "eBMD", "FEV1FVC", "Fibroblastic_Disorders", "GGT", "Glaucoma_Combined", "Inguinal_Hernia", "LOY", "MCP")
)
domain_rows <- rbindlist(lapply(names(domain_map), function(d) {
  data.table(domain = d, trait = unique(domain_map[[d]]))
}), fill = TRUE)
domain_rows <- unique(domain_rows)
fwrite(domain_rows, file.path(out, "config", "ukbb_trait_domains.tsv"), sep = "\t")

domain_lookup <- setNames(domain_rows$domain, domain_rows$trait)
uk[, domain := unname(domain_lookup[as.character(trait)])]
uk[is.na(domain), domain := "other"]

all_traits <- sort(unique(uk$trait))
fwrite(data.table(trait = all_traits, n_rows = vapply(all_traits, function(t) sum(uk$trait == t), integer(1))),
       file.path(out, "config", "ukbb_traits_observed.tsv"), sep = "\t")

bootstrap_mean <- function(x, n_boot, seed) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(numeric())
  set.seed(seed)
  replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
}

summarize_trait_set <- function(x, analysis, trait_set, bootstrap_unit = "trait") {
  x <- x[trait %in% trait_set]
  x[, analysis := analysis]
  ans <- x[, {
    valid <- is.finite(score)
    vals <- score[valid]
    traits_valid <- unique(trait[valid])
    boot <- bootstrap_mean(vals, n_boot, seed + .GRP)
    list(
      n_selected = length(trait_set), n_valid = length(vals),
      mean_log2_1p = if (length(vals)) mean(vals) else NA_real_,
      median_log2_1p = if (length(vals)) median(vals) else NA_real_,
      geomean_1p_minus1 = if (length(vals)) 2^mean(vals) - 1 else NA_real_,
      ci_low_log2_1p = if (length(boot)) quantile(boot, 0.025, names = FALSE) else NA_real_,
      ci_high_log2_1p = if (length(boot)) quantile(boot, 0.975, names = FALSE) else NA_real_,
      ci_low_enrichment = if (length(boot)) 2^quantile(boot, 0.025, names = FALSE) - 1 else NA_real_,
      ci_high_enrichment = if (length(boot)) 2^quantile(boot, 0.975, names = FALSE) - 1 else NA_real_,
      bootstrap_replicates = if (length(boot)) n_boot else NA_integer_,
      bootstrap_unit = bootstrap_unit,
      traits_valid = paste(traits_valid, collapse = ";")
    )
  }, by = .(analysis, view, dataset, celltype, dataset_celltype, matched_panel, method_label)]
  ans[]
}

primary_parts <- lapply(names(panels), function(p) {
  x <- uk[matched_panel == p]
  summarize_trait_set(x, "primary_biology_matched", panels[[p]], "trait")
})
primary <- rbindlist(primary_parts, fill = TRUE)

# Supplementary all-94: same log-scale metric, but all traits are used.
supplementary <- summarize_trait_set(uk, "supplementary_all94", all_traits, "trait")

# Secondary domain-balanced result. First average within each domain, then
# average the domain means. Bootstrap samples domains, not individual traits,
# so domains receive equal weight.
domain_trait_summary <- uk[is.finite(score), .(
  domain_score = mean(score), n_domain_traits = .N,
  domain_traits = paste(sort(unique(trait)), collapse = ";")
), by = .(view, dataset, celltype, dataset_celltype, matched_panel, method_label, domain)]

domain_balanced <- domain_trait_summary[, {
  vals <- domain_score[is.finite(domain_score)]
  boot <- if (length(vals) >= 2L) {
    set.seed(seed + .GRP)
    replicate(n_boot, mean(sample(vals, length(vals), replace = TRUE)))
  } else numeric()
  list(
    n_domains = length(vals), n_traits = sum(n_domain_traits),
    mean_log2_1p = if (length(vals)) mean(vals) else NA_real_,
    median_log2_1p = if (length(vals)) median(vals) else NA_real_,
    geomean_1p_minus1 = if (length(vals)) 2^mean(vals) - 1 else NA_real_,
    ci_low_log2_1p = if (length(boot)) quantile(boot, .025, names = FALSE) else NA_real_,
    ci_high_log2_1p = if (length(boot)) quantile(boot, .975, names = FALSE) else NA_real_,
    ci_low_enrichment = if (length(boot)) 2^quantile(boot, .025, names = FALSE) - 1 else NA_real_,
    ci_high_enrichment = if (length(boot)) 2^quantile(boot, .975, names = FALSE) - 1 else NA_real_,
    bootstrap_replicates = if (length(boot)) n_boot else NA_integer_,
    bootstrap_unit = "domain",
    domains = paste(sort(unique(domain[is.finite(domain_score)])), collapse = ";")
  )
}, by = .(view, dataset, celltype, dataset_celltype, matched_panel, method_label)]
domain_balanced[, analysis := "secondary_domain_balanced"]

write_table <- function(x, path) {
  y <- copy(x)
  for (nm in intersect(c("method_label"), names(y))) y[, (nm) := as.character(get(nm))]
  fwrite(y, path, sep = "\t", na = "NA")
}

write_table(primary, file.path(out, "tables", "primary_biology_matched.tsv"))
write_table(domain_balanced, file.path(out, "tables", "secondary_domain_balanced.tsv"))
write_table(supplementary, file.path(out, "tables", "supplementary_all94.tsv"))
write_table(domain_trait_summary, file.path(out, "tables", "secondary_domain_level.tsv"))

# Paired method comparisons on the same traits. Positive delta means method_b
# has larger log2(1 + enrichment) than method_a for the same trait.
pair_rows <- list(); pk <- 0L
pair_methods <- combn(methods, 2L, simplify = FALSE)
pair_input <- uk[, .(analysis_group = matched_panel, view, dataset, celltype,
                     dataset_celltype, trait, method_label = as.character(method_label), score)]
for (pp in pair_methods) {
  a <- pair_input[method_label == pp[[1L]], .(view, dataset, celltype, dataset_celltype,
                                              trait, score_a = score)]
  b <- pair_input[method_label == pp[[2L]], .(view, dataset, celltype, dataset_celltype,
                                              trait, score_b = score)]
  z <- merge(a, b, by = c("view", "dataset", "celltype", "dataset_celltype", "trait"), all = FALSE)
  z[, delta := score_b - score_a]
  z <- z[is.finite(delta)]
  if (!nrow(z)) next
  z[, matched_panel := mapply(panel_map, dataset, celltype)]
  groups <- unique(z[, .(view, dataset, celltype, dataset_celltype, matched_panel)])
  for (ii in seq_len(nrow(groups))) {
    ag <- groups[ii]
    q0 <- z[view == ag$view & dataset == ag$dataset & celltype == ag$celltype]
    for (analysis_set in c("primary_biology_matched", "supplementary_all94")) {
      q <- if (analysis_set == "primary_biology_matched")
        q0[trait %in% panels[[ag$matched_panel]]] else q0
      vals <- q$delta
      set.seed(seed + pk + 1L)
      boot <- if (length(vals) >= 2L) replicate(n_boot, mean(sample(vals, length(vals), replace = TRUE))) else numeric()
      pval <- if (length(vals) >= 2L && any(vals != 0)) suppressWarnings(wilcox.test(vals, mu = 0, exact = FALSE)$p.value) else NA_real_
      pk <- pk + 1L
      pair_rows[[pk]] <- data.table(
        analysis = analysis_set, view = ag$view, dataset = ag$dataset,
        celltype = ag$celltype, dataset_celltype = ag$dataset_celltype,
        matched_panel = ag$matched_panel, method_a = pp[[1L]], method_b = pp[[2L]],
        n_traits = length(vals), mean_delta_log2_1p = if (length(vals)) mean(vals) else NA_real_,
        median_delta_log2_1p = if (length(vals)) median(vals) else NA_real_,
        win_fraction_method_b = if (length(vals)) mean(vals > 0) else NA_real_,
        ci_low_delta = if (length(boot)) quantile(boot, .025, names = FALSE) else NA_real_,
        ci_high_delta = if (length(boot)) quantile(boot, .975, names = FALSE) else NA_real_,
        paired_wilcoxon_p = pval, bootstrap_replicates = if (length(boot)) n_boot else NA_integer_,
        bootstrap_unit = "trait"
      )
    }
  }
}
paired <- rbindlist(pair_rows, fill = TRUE)
write_table(paired, file.path(out, "tables", "paired_method_comparisons.tsv"))

write_table(uk[, .(source, annotation_id, dataset, celltype, view, matched_panel, domain,
                   method_label, trait, trait_label, enrichment, score, n_causal_total,
                   n_causal_in_annotation, n_common_locus_total, n_common_in_annotation,
                   fisher_p, source_file)],
            file.path(out, "tables", "per_trait_input_used.tsv"))

theme_panel <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))

plot_summary <- function(x, file, title, subtitle) {
  if (!nrow(x)) return(invisible(NULL))
  x <- copy(x)
  x[, method_label := factor(as.character(method_label), levels = methods)]
  x[, row_label := ifelse(view == "pooled", dataset, dataset_celltype)]
  x[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  p <- ggplot(x, aes(mean_log2_1p, row_label, color = method_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = .35) +
    geom_errorbar(aes(xmin = ci_low_log2_1p, xmax = ci_high_log2_1p),
                  orientation = "y", height = .18, position = position_dodge(width = .65),
                  linewidth = .45, na.rm = TRUE) +
    geom_point(position = position_dodge(width = .65), size = 1.9, na.rm = TRUE) +
    facet_wrap(~matched_panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = method_colors, drop = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = "Mean log2(1 + enrichment)", y = NULL, color = "Method") + theme_panel
  h <- max(5.5, .28 * length(unique(x$row_label)) + 2.2)
  ggsave(file, p, width = 11.5, height = h, units = "in", device = cairo_pdf)
}

for (vw in c("pooled", "celltype")) {
  plot_summary(primary[view == vw], file.path(out, "figures", paste0("primary_", vw, ".pdf")),
               "Biology-matched UKBB SCENT enrichment", paste0(vw, "; point = mean log2(1 + enrichment), whisker = 1,000-trait bootstrap CI"))
  plot_summary(domain_balanced[view == vw], file.path(out, "figures", paste0("secondary_domain_balanced_", vw, ".pdf")),
               "Domain-balanced UKBB SCENT enrichment", paste0(vw, "; each predefined trait domain has equal weight; whisker = domain bootstrap CI"))
  plot_summary(supplementary[view == vw], file.path(out, "figures", paste0("supplementary_all94_", vw, ".pdf")),
               "Supplementary all-94-trait UKBB SCENT enrichment", paste0(vw, "; mean log2(1 + enrichment), whisker = 1,000-trait bootstrap CI"))
}

plot_paired <- function(vw, file) {
  z <- paired[view == vw & analysis == "primary_biology_matched"]
  if (!nrow(z)) return(invisible(NULL))
  z[, comparison := paste(method_b, "vs", method_a)]
  z[, comparison := factor(comparison,
                           levels = unique(paste(method_b, "vs", method_a)))]
  z[, row_label := ifelse(view == "pooled", dataset, dataset_celltype)]
  z[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  p <- ggplot(z, aes(comparison, row_label, fill = mean_delta_log2_1p)) +
    geom_tile(color = "white", linewidth = .25) +
    geom_text(aes(label = sprintf("%+.2f", mean_delta_log2_1p)), size = 2.5) +
    facet_wrap(~matched_panel, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = "Mean paired delta\nlog2(1 + E)") +
    labs(title = "Pairwise method comparison on the same UKBB traits",
         subtitle = paste0(vw, "; positive values favor the method listed first"),
         x = NULL, y = NULL) + theme_panel
  h <- max(5.5, .28 * length(unique(z$row_label)) + 2.2)
  ggsave(file, p, width = 10.5, height = h, units = "in", device = cairo_pdf)
}
plot_paired("pooled", file.path(out, "figures", "paired_primary_pooled.pdf"))
plot_paired("celltype", file.path(out, "figures", "paired_primary_celltype.pdf"))

manifest <- data.table(
  output = c("primary", "secondary", "supplementary"),
  definition = c(
    "Dataset/cell-type matched trait panel; arithmetic mean of log2(1 + enrichment); bootstrap unit = trait",
    "All observed UKBB traits, averaged within predefined domains and then equally across domains; bootstrap unit = domain",
    "All observed UKBB traits; arithmetic mean of log2(1 + enrichment); bootstrap unit = trait"
  ),
  n_boot = c(n_boot, n_boot, n_boot), seed = c(seed, seed, seed),
  methods = paste(methods, collapse = "; ")
)
fwrite(manifest, file.path(out, "config", "analysis_manifest.tsv"), sep = "\t")

message("Wrote UKBB trait-panel benchmark to ", out)
message("Rows: primary=", nrow(primary), "; domain-balanced=", nrow(domain_balanced),
        "; all94=", nrow(supplementary), "; paired=", nrow(paired))
