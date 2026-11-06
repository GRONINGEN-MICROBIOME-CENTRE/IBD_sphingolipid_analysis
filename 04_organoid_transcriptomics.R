source("public_code/00_plot_helpers.R")
library(broom)

# expression: sample, gene, expression (logCPM), cytokines, timepoint,
# donor_line, valeric_acid; tested_genes: selected 67-gene testing set.
# membership: gene, module, order (51 main-panel genes).
make_organoid_figures <- function(expression, tested_genes, membership,
                                  out_dir = "figures/organoids") {
  stopifnot(length(unique(tested_genes)) == 67, nrow(membership) == 51,
            !anyDuplicated(membership$gene), all(membership$gene %in% tested_genes))
  d <- expression |> filter(valeric_acid == 0, gene %in% tested_genes) |>
    mutate(cytokines = factor(cytokines, levels = c("Control", "Cytokines")),
           timepoint = factor(timepoint), donor_line = factor(donor_line))
  fit_effect <- function(d) {
    fit <- lm(value ~ cytokines + timepoint + donor_line, data = d)
    tidy(fit) |> filter(term == "cytokinesCytokines")
  }
  genes <- d |> mutate(value = expression) |> group_by(gene) |>
    group_modify(~ fit_effect(.x)) |> ungroup() |> mutate(fdr = p.adjust(p.value, "BH"))
  heat <- genes |> inner_join(membership, by = "gene") |>
    mutate(gene = factor(gene, levels = rev(membership$gene[order(membership$order)])),
           module = factor(module, levels = module_levels), star = stars(fdr))
  hp <- ggplot(heat, aes("Cytokines vs control", gene, fill = estimate)) +
    geom_tile(colour = "white") + geom_text(aes(label = star), size = 3) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    labs(x = NULL, y = NULL, fill = "Adjusted\nlogCPM difference") + theme_figure() +
    theme(axis.text.y = element_text(size = 7), strip.text.y = element_text(angle = 0))
  scores <- d |> inner_join(membership, by = "gene") |> group_by(gene) |>
    mutate(z = safe_z(expression)) |> ungroup() |>
    group_by(sample, module, cytokines, timepoint, donor_line) |>
    summarise(value = mean(z), .groups = "drop")
  module_tests <- scores |> group_by(module) |> group_modify(~ fit_effect(.x)) |>
    ungroup() |> mutate(fdr = p.adjust(p.value, "BH"))
  # Display adjustment follows the source module panels: remove timepoint and line effects.
  scores <- scores |> group_by(module) |> group_modify(function(d, key) {
    fit <- lm(value ~ timepoint + donor_line, data = d)
    mutate(d, adjusted_score = resid(fit) + mean(value))
  }) |> ungroup()
  selected <- module_levels[c(1,2,3,5)]
  panels <- lapply(selected, function(m) {
    s <- filter(scores, module == m); r <- filter(module_tests, module == m)
    ggplot(s, aes(cytokines, adjusted_score, fill = cytokines)) +
      geom_violin(trim = FALSE, colour = NA, alpha = .65) +
      geom_boxplot(width = .2, outlier.shape = NA) +
      geom_point(position = position_jitter(width = .07, seed = 7), size = 1.5,
                 colour = "grey35") +
      annotate("text", x = 1.5, y = max(s$adjusted_score) + .15 * diff(range(s$adjusted_score)),
               label = p_label(r$fdr), size = 3) +
      scale_fill_manual(values = c(Control = "#6494C1", Cytokines = "#D15E60"), guide = "none") +
      labs(x = NULL, y = "Adjusted module score", title = m) + theme_figure()
  })
  save_figure(hp, out_dir, "Figure5A", 6, 11)
  save_figure(wrap_plots(panels, ncol = 2), out_dir, "Figure5B", 10, 8)
  write_csv(genes, file.path(out_dir, "Figure5_gene_effects.csv"))
  write_csv(module_tests, file.path(out_dir, "Figure5_module_effects.csv"))
}
