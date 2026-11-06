library(tidyverse)
library(patchwork)

module_levels <- c("De novo ceramide synthesis", "Complex sphingolipid synthesis",
                   "Degradation and salvage", "S1P metabolism",
                   "S1P receptors and signaling", "Transport and trafficking",
                   "Regulation and auxiliary branches")
class_levels <- c("Ceramide", "Sphingomyelin",
                  "N-acyl/glycosyl/lactosyl sphingolipid", "Sphingoid base",
                  "Other sphingolipid-like")
location_levels <- c("Control", "Colonic IBD", "Ileal IBD")
status_levels <- c("Control", "IBD non-inflamed", "IBD inflamed")

safe_z <- function(x) {
  if (sum(is.finite(x)) < 2 || isTRUE(sd(x, na.rm = TRUE) == 0)) return(rep(0, length(x)))
  as.numeric(scale(x))
}
half_min_log <- function(x, base = 2) {
  positive <- x[is.finite(x) & x > 0]
  pseudo <- if (length(positive)) min(positive) / 2 else 1
  log(x + pseudo, base = base)
}
p_label <- function(p, prefix = "FDR") {
  ifelse(is.na(p), paste0(prefix, "=NA"),
         ifelse(p < .001, paste0(prefix, "<0.001"),
                paste0(prefix, "=", formatC(p, format = "f", digits = 3))))
}
stars <- function(p) case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
theme_figure <- function() {
  theme_classic(base_size = 10) +
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text = element_text(face = "bold"))
}
save_figure <- function(plot, out_dir, name, width, height) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height,
         limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), plot, width = width, height = height,
         dpi = 300, limitsize = FALSE)
}
location_group <- function(diagnosis, disease_location) {
  case_when(diagnosis == "Control" ~ "Control",
            diagnosis %in% c("CD", "UC", "IBDU", "IBD") & disease_location == "colon" ~ "Colonic IBD",
            diagnosis %in% c("CD", "UC", "IBDU", "IBD") & disease_location == "ileum" ~ "Ileal IBD",
            TRUE ~ NA_character_)
}
# Input: feature, value, group. Contrast estimates are group2 minus group1.
rank_statistics <- function(dat, contrasts) {
  dat <- filter(dat, is.finite(value), !is.na(group))
  overall <- dat |> group_by(feature) |> summarise(
    p_value = if (n_distinct(group) > 1) kruskal.test(value ~ group)$p.value else NA_real_,
    .groups = "drop") |> mutate(fdr = p.adjust(p_value, "BH"))
  pairwise <- dat |> group_by(feature) |> group_modify(function(d, key) {
    map_dfr(seq_len(nrow(contrasts)), function(i) {
      x <- d$value[d$group == contrasts$group1[i]]
      y <- d$value[d$group == contrasts$group2[i]]
      tibble(contrast = contrasts$contrast[i], group1 = contrasts$group1[i],
             group2 = contrasts$group2[i],
             estimate = median(y) - median(x),
             p_value = if (length(x) >= 2 && length(y) >= 2)
               suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value) else NA_real_)
    })
  }) |> ungroup() |> group_by(contrast) |>
    mutate(fdr = p.adjust(p_value, "BH")) |> ungroup() |>
    mutate(fdr_global = p.adjust(p_value, "BH"))
  list(overall = overall, pairwise = pairwise)
}
all_contrasts <- function(levels) {
  pairs <- combn(levels, 2)
  tibble(group1 = pairs[1, ], group2 = pairs[2, ],
         contrast = paste(group1, group2, sep = " vs "))
}
# Statistics supplied here must already have been corrected across all features.
distribution_plot <- function(dat, results, group_levels, y_label, raw_p = FALSE) {
  map_plots <- lapply(unique(dat$feature), function(f) {
    d <- filter(dat, feature == f, is.finite(value)) |>
      mutate(group = factor(group, levels = group_levels))
    r <- filter(results, feature == f)
    span <- diff(range(d$value)); if (!is.finite(span) || span == 0) span <- 1
    r <- r |> mutate(x1 = match(group1, group_levels), x2 = match(group2, group_levels),
                     y = max(d$value) + seq_len(n()) * .13 * span,
                     label = p_label(if (raw_p) p_value else fdr, if (raw_p) "P" else "FDR"))
    ggplot(d, aes(group, value, fill = group)) +
      geom_violin(trim = FALSE, alpha = .65, colour = NA) +
      geom_boxplot(width = .2, outlier.shape = NA) +
      geom_point(position = position_jitter(width = .08, seed = 7), size = .7,
                 alpha = .4, colour = "grey35") +
      geom_segment(data = r, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE) +
      geom_text(data = r, aes(x = (x1 + x2)/2, y = y, label = label),
                inherit.aes = FALSE, vjust = -.3, size = 2.7) +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
      labs(x = NULL, y = y_label, title = f) + theme_figure()
  })
  wrap_plots(map_plots, ncol = 2)
}
write_statistics <- function(results, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(results$overall, file.path(out_dir, paste0(prefix, "_overall.csv")))
  write_csv(results$pairwise, file.path(out_dir, paste0(prefix, "_contrasts.csv")))
}
