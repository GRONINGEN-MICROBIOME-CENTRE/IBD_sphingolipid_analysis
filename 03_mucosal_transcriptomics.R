source("public_code/00_plot_helpers.R")
library(lmerTest)
library(emmeans)

# expression: sample, participant, site, group, gene, expression (ComBat logCPM).
# Optional covariates: age, sex, bmi. membership: gene, module, order (52 genes).
fit_mucosal_feature <- function(d, adjusted = FALSE) {
  d$group <- factor(d$group, levels = status_levels)
  d$participant <- factor(d$participant)
  if (adjusted) d$sex <- factor(d$sex)
  formula <- if (adjusted) value ~ group + age + sex + bmi + (1 | participant) else
    value ~ group + (1 | participant)
  fit <- lmer(formula, data = d, REML = TRUE, na.action = na.omit)
  emm <- emmeans(fit, ~ group, lmer.df = "satterthwaite")
  comparisons <- contrast(emm, list(
    noninflamed_minus_control = c(-1,1,0), inflamed_minus_control = c(-1,0,1),
    inflamed_minus_noninflamed = c(0,-1,1)), adjust = "none") |>
    summary(infer = c(TRUE,TRUE)) |> as_tibble() |>
    mutate(group1 = c(status_levels[1], status_levels[1], status_levels[2]),
           group2 = c(status_levels[2], status_levels[3], status_levels[3]))
  list(means = as_tibble(summary(emm, infer = c(TRUE,FALSE))), contrasts = comparisons)
}
mucosal_models <- function(d, adjusted = FALSE) {
  keys <- distinct(d, site, feature)
  models <- map(seq_len(nrow(keys)), function(i) {
    fit_mucosal_feature(filter(d, site == keys$site[i], feature == keys$feature[i]), adjusted)
  })
  means <- map_dfr(seq_along(models), function(i) bind_cols(keys[i, ], models[[i]]$means))
  tests <- map_dfr(seq_along(models), function(i) bind_cols(keys[i, ], models[[i]]$contrasts)) |>
    group_by(site, contrast) |> mutate(fdr = p.adjust(p.value, "BH")) |> ungroup()
  list(means = means, tests = tests)
}
make_mucosal_figures <- function(expression, membership, out_dir = "figures/mucosal",
                                 sensitivity = FALSE) {
  stopifnot(nrow(membership) == 52, !anyDuplicated(membership$gene))
  d <- expression |> inner_join(membership, by = "gene") |>
    filter(is.finite(expression), !is.na(site), !is.na(group), !is.na(participant)) |>
    mutate(site = factor(site, levels = c("Colon", "Ileum")),
           group = factor(group, levels = status_levels))
  genes <- d |> mutate(feature = gene, value = expression)
  modules <- d |> group_by(gene) |> mutate(z = safe_z(expression)) |> ungroup() |>
    group_by(sample, participant, site, group, module) |>
    summarise(value = mean(z), .groups = "drop") |> mutate(feature = module)
  gene_fit <- mucosal_models(genes)
  module_fit <- mucosal_models(modules)
  heat <- gene_fit$means |> rename(gene = feature) |> inner_join(membership, by = "gene") |>
    group_by(gene) |> mutate(z = safe_z(emmean)) |> ungroup() |>
    mutate(gene = factor(gene, levels = rev(membership$gene[order(membership$order)])),
           module = factor(module, levels = module_levels))
  hp <- ggplot(heat, aes(group, gene, fill = z)) + geom_tile(colour = "white", linewidth = .24) +
    facet_grid(module ~ site, scales = "free", space = "free") +
    scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#B83232",
                         limits = c(-2.5,2.5), oob = scales::squish) +
    labs(x = NULL, y = NULL, fill = "Gene-wise z-score") + theme_figure() +
    theme(axis.text.y = element_text(size = 6), strip.text.y = element_text(angle = 0))
  save_figure(hp, out_dir, "Figure4A", 9, 12)
  selected_modules <- module_levels[c(1,2,3,5)]
  panels <- lapply(selected_modules, function(m) {
    # BH values are calculated across all seven modules before panel selection.
    plots <- lapply(c("Colon", "Ileum"), function(s) {
      distribution_plot(filter(modules, feature == m, site == s),
                        filter(module_fit$tests, feature == m, site == s) |>
                          rename(p_value = p.value), status_levels, "Mean gene z-score") +
        plot_annotation(title = paste(m, s, sep = ": "))
    })
    wrap_plots(plots, ncol = 2)
  })
  save_figure(wrap_plots(panels, ncol = 1), out_dir, "Figure4B", 13, 15)
  write_csv(gene_fit$tests, file.path(out_dir, "Figure4_gene_contrasts.csv"))
  write_csv(module_fit$tests, file.path(out_dir, "Figure4_module_contrasts.csv"))
  if (sensitivity) {
    covariates <- d |> distinct(sample, age, sex, bmi) |>
      filter(is.finite(age), !is.na(sex), is.finite(bmi))
    for (layer in c("genes", "modules")) {
      values <- if (layer == "genes") genes else modules
      complete <- values |> select(-any_of(c("age", "sex", "bmi"))) |>
        inner_join(covariates, by = "sample")
      baseline <- mucosal_models(complete)$tests
      adjusted <- mucosal_models(complete, adjusted = TRUE)$tests
      comparison <- adjusted |> select(site, feature, contrast, adjusted = estimate, adjusted_fdr = fdr) |>
        inner_join(baseline |> select(site, feature, contrast, baseline = estimate, baseline_fdr = fdr),
                   by = c("site", "feature", "contrast"))
      plot <- ggplot(comparison, aes(baseline, adjusted)) + geom_abline(slope = 1, intercept = 0,
        colour = "grey60", linetype = 2) + geom_point(size = 1.5) +
        facet_grid(contrast ~ site, scales = "free") + theme_figure() +
        labs(x = "Baseline effect (complete cases)", y = "Age, sex and BMI adjusted effect")
      save_figure(plot, out_dir, paste0("Figure4_sensitivity_", layer), 10, 9)
      write_csv(comparison, file.path(out_dir, paste0("Figure4_sensitivity_", layer, ".csv")))
    }
  }
}
