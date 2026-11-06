source("public_code/00_plot_helpers.R")

# mag_abundance: sample, mag, relative_abundance; detected MAGs only.
# retained_hits: mag, target_protein, bitscore; metadata: sample, diagnosis,
# disease_location, fcp. Input hits have passed sequence-search filtering.
make_mag_figures <- function(mag_abundance, retained_hits, metadata,
                             out_dir = "figures/metagenomics") {
  hits <- retained_hits |> group_by(mag, target_protein) |>
    slice_max(bitscore, n = 1, with_ties = FALSE) |> ungroup()
  positive_mags <- unique(hits$mag)
  metrics <- mag_abundance |> mutate(hit = mag %in% positive_mags) |>
    group_by(sample) |> summarise(
      relative_abundance = sum(relative_abundance[hit], na.rm = TRUE),
      proportion = n_distinct(mag[hit]) / n_distinct(mag), .groups = "drop") |>
    mutate(relative_abundance = log1p(relative_abundance)) |>
    inner_join(metadata, by = "sample") |>
    mutate(group3 = location_group(diagnosis, disease_location),
           fcp_group = case_when(is.finite(fcp) & fcp > 200 ~ "FCP > 200",
                                 is.finite(fcp) & fcp <= 200 ~ "FCP <= 200"),
           group5 = case_when(group3 == "Control" ~ "Control",
                              !is.na(group3) & !is.na(fcp_group) ~ paste(group3, fcp_group, sep = "\n"))) |>
    pivot_longer(c(relative_abundance, proportion), names_to = "feature", values_to = "value")
  for (panel in c("A", "B")) {
    levels <- if (panel == "A") location_levels else c("Control",
      "Colonic IBD\nFCP <= 200", "Colonic IBD\nFCP > 200",
      "Ileal IBD\nFCP <= 200", "Ileal IBD\nFCP > 200")
    column <- if (panel == "A") "group3" else "group5"
    d <- metrics |> mutate(group = factor(.data[[column]], levels = levels)) |>
      filter(!is.na(group), is.finite(value))
    tests <- rank_statistics(d, all_contrasts(levels))
    plot <- distribution_plot(d, tests$pairwise, levels, "SL-hit MAG metric", raw_p = TRUE)
    save_figure(plot, out_dir, paste0("Figure3", panel), 12, if (panel == "A") 5 else 8)
    write_statistics(tests, out_dir, paste0("Figure3", panel))
  }
}

# features: sample, uniref, abundance; mapping: uniref, enzyme;
# metabolites: sample, feature, abundance, module (paired IBD profiles).
# enzyme_order: enzyme names ordered by the separate IBD/control mean ratio.
make_enzyme_heatmap <- function(features, mapping, metabolites, enzyme_order,
                                out_dir = "figures/metagenomics") {
  mapped <- features |> inner_join(distinct(mapping), by = "uniref")
  positive <- mapped$abundance[is.finite(mapped$abundance) & mapped$abundance > 0]
  stopifnot(length(positive) > 0)
  pseudo <- min(positive) / 2
  paired <- intersect(unique(features$sample), unique(metabolites$sample))
  enzymes <- mapped |> filter(sample %in% paired) |>
    mutate(abundance = replace_na(abundance, 0)) |>
    group_by(sample, enzyme) |> summarise(abundance = sum(abundance), .groups = "drop") |>
    complete(sample = paired, enzyme = unique(mapping$enzyme), fill = list(abundance = 0)) |>
    mutate(log_enzyme = log(abundance + pseudo))
  lipids <- metabolites |> filter(sample %in% paired) |> group_by(feature) |>
    mutate(log_metabolite = half_min_log(abundance, base = exp(1))) |> ungroup()
  joined <- inner_join(lipids, enzymes, by = "sample", relationship = "many-to-many")
  results <- joined |> group_by(feature, module, enzyme) |> group_modify(function(d, key) {
    d <- filter(d, is.finite(log_metabolite), is.finite(log_enzyme))
    if (nrow(d) < 10 || sd(d$log_metabolite) == 0 || sd(d$log_enzyme) == 0)
      return(tibble(rho = NA_real_, p_value = NA_real_))
    test <- suppressWarnings(cor.test(d$log_metabolite, d$log_enzyme,
                                     method = "spearman", exact = FALSE))
    tibble(rho = unname(test$estimate), p_value = test$p.value)
  }) |> ungroup() |> mutate(fdr = p.adjust(p_value, "BH"))
  display <- results |> filter(is.finite(rho), enzyme %in% enzyme_order) |>
    arrange(fdr, p_value) |> slice_head(n = 80) |>
    mutate(enzyme = factor(enzyme, levels = enzyme_order),
           module = factor(module, levels = class_levels), star = stars(fdr))
  plot <- ggplot(display, aes(enzyme, feature, fill = rho)) +
    geom_tile(colour = "white") + geom_text(aes(label = star), size = 3) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-.5,.5), oob = scales::squish) +
    labs(x = NULL, y = NULL, fill = "Spearman rho") + theme_figure() +
    theme(axis.text.y = element_text(size = 6), strip.text.y = element_text(angle = 0))
  save_figure(plot, out_dir, "Figure3C", 10, 11)
  write_csv(results, file.path(out_dir, "Figure3C_correlations.csv"))
}
