source("public_code/00_plot_helpers.R")

# abundance: sample, feature, abundance; metadata: sample, diagnosis,
# disease_location, fcp; membership: feature, module.
# membership defines the curated 46-metabolite panel and five classes.
make_metabolomics_figures <- function(abundance, metadata, membership,
                                      out_dir = "figures/metabolomics") {
  stopifnot(!anyDuplicated(membership$feature), !anyDuplicated(metadata$sample))
  d <- abundance |> inner_join(membership, by = "feature") |>
    inner_join(metadata, by = "sample") |>
    mutate(group3 = location_group(diagnosis, disease_location),
           fcp_group = case_when(is.finite(fcp) & fcp > 200 ~ "FCP > 200",
                                 is.finite(fcp) & fcp <= 200 ~ "FCP <= 200"),
           group6 = if_else(!is.na(group3) & !is.na(fcp_group),
                            paste(group3, fcp_group, sep = "\n"), NA_character_)) |>
    filter(!is.na(group3), is.finite(abundance))
  sums <- bind_rows(
    d |> group_by(sample, group3, group6, feature = module) |>
      summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop"),
    d |> group_by(sample, group3, group6) |>
      summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") |>
      mutate(feature = "Total sphingolipids")) |>
    group_by(feature) |> mutate(value = half_min_log(abundance)) |> ungroup()
  for (fig in 1:2) {
    levels <- if (fig == 1) location_levels else
      as.vector(rbind(paste(location_levels, "FCP <= 200", sep = "\n"),
                      paste(location_levels, "FCP > 200", sep = "\n")))
    column <- if (fig == 1) "group3" else "group6"
    comparisons <- if (fig == 1) all_contrasts(levels) else tibble(
      group1 = levels[c(1,3,5)], group2 = levels[c(2,4,6)],
      contrast = paste(location_levels, "high vs low FCP"))
    metabolites <- d |> mutate(group = factor(.data[[column]], levels = levels),
                                value = log2(abundance + 1)) |> filter(!is.na(group))
    tests <- rank_statistics(metabolites, comparisons)
    selected <- tests$pairwise |> filter(!is.na(fdr), fdr < .05) |> pull(feature)
    if (fig == 1) selected <- union(selected, filter(tests$overall, fdr < .05)$feature)
    medians <- metabolites |> group_by(feature, module, group) |>
      summarise(median_value = median(value), .groups = "drop") |>
      group_by(feature) |> mutate(z = safe_z(median_value)) |> ungroup()
    row_order <- membership |> filter(feature %in% selected) |>
      mutate(module = factor(module, levels = class_levels)) |> arrange(module, feature)
    heat <- medians |> filter(feature %in% selected) |>
      mutate(feature = factor(feature, levels = rev(row_order$feature)),
             module = factor(module, levels = class_levels))
    hp <- ggplot(heat, aes(group, feature, fill = z)) +
      geom_tile(colour = "white", linewidth = .25) +
      facet_grid(module ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B43C3C",
                           limits = c(-2,2), oob = scales::squish) +
      labs(x = NULL, y = NULL, fill = "Median\nrow z-score") + theme_figure() +
      theme(axis.text.y = element_text(size = 6), strip.text.y = element_text(angle = 0))
    if (nrow(heat) == 0) hp <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No metabolites pass FDR < 0.05") + theme_void()
    modules <- sums |> mutate(group = factor(.data[[column]], levels = levels)) |>
      filter(!is.na(group))
    module_tests <- rank_statistics(modules, comparisons)
    mp <- distribution_plot(modules, module_tests$pairwise, levels, "log2 summed abundance")
    save_figure(hp, out_dir, paste0("Figure", fig, "A"), if (fig == 1) 9 else 11, 10)
    save_figure(mp, out_dir, paste0("Figure", fig, "B"), 11, 12)
    write_statistics(tests, out_dir, paste0("Figure", fig, "_metabolites"))
    write_statistics(module_tests, out_dir, paste0("Figure", fig, "_modules"))
  }
}
