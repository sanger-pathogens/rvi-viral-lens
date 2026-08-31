#!/usr/bin/env Rscript
# Heatmap figure: relative abundance before/after SCRuB decontamination, and the change.
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
})

CONTROL_COLOR <- "#c0392b"
MIN_READ_CHANGE <- 20  # below this raw read-count change, treat the cell as unchanged (not colored)

option_list <- list(
  make_option("--bracken-tsv", type = "character", dest = "bracken_tsv",
              help = "Reformatted samples x species Bracken TSV (SCRuB's input matrix)"),
  make_option("--scrub-rds", type = "character", dest = "scrub_rds",
              help = "scrub_result.rds produced by the SCRUB process"),
  make_option("--top-n", type = "integer", default = 40, dest = "top_n",
              help = "Number of most abundant species to plot [default %default]"),
  make_option("--out", type = "character", default = "decontam_heatmap.png", dest = "out",
              help = "Output PNG path [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

## ---- Load data ----------------------------------------------------------
counts_before_raw <- as.matrix(read.delim(opt$bracken_tsv, row.names = 1, check.names = FALSE))
scrub <- readRDS(opt$scrub_rds)
counts_after_raw <- scrub$decontaminated_samples

## SCRuB drops the control/blank samples it uses internally when it returns the
## decontaminated matrix. Keep them for the "before" panel (flagged in red) but
## they have no "after"/"change" data.
control_samples <- setdiff(rownames(counts_before_raw), rownames(counts_after_raw))
if (length(control_samples) > 0) {
  message("Samples absent from decontaminated output (SCRuB control/blank samples): ",
          paste(control_samples, collapse = ", "))
}

## Order ALL samples numerically by the "#N" suffix (e.g. 50376_2#1 ... 50376_2#96)
## and reuse this order + factor levels across all three panels so rows line up.
sample_levels <- rownames(counts_before_raw)[order(as.numeric(sub(".*#", "", rownames(counts_before_raw))))]
is_control <- setNames(sample_levels %in% control_samples, sample_levels)
after_order <- sample_levels[!is_control[sample_levels]]

## SCRuB re-read the species names through R's make.names() (hyphens/slashes -> dots),
## but column order is preserved, so align positionally and keep the original labels.
species <- colnames(counts_before_raw)
stopifnot(identical(make.names(species), colnames(counts_after_raw)))
colnames(counts_after_raw) <- species

counts_before <- counts_before_raw[sample_levels, species]
counts_after  <- counts_after_raw[after_order, species]

## ---- Relative abundance ---------------------------------------------------
row_norm <- function(m) {
  rs <- rowSums(m)
  rs[rs == 0] <- 1  # avoid 0/0; rows with no reads stay all-zero
  m / rs
}

relabund_before <- row_norm(counts_before)                          # all samples, incl. controls
relabund_after  <- row_norm(counts_after)                           # decontaminated samples only
change          <- relabund_after - row_norm(counts_before[after_order, species])

## Raw read-count change is noisy at low counts; mask out cells where SCRuB moved fewer
## than MIN_READ_CHANGE reads so the change panel only colors meaningful shifts.
read_count_change <- counts_after - counts_before[after_order, species]
change[abs(read_count_change) < MIN_READ_CHANGE] <- 0

## ---- Pick the top N most abundant species (by mean relative abundance pre-decontam) ----
top_species <- names(sort(colMeans(relabund_before), decreasing = TRUE))[seq_len(min(opt$top_n, ncol(relabund_before)))]

to_long <- function(m, value_name) {
  m[, top_species, drop = FALSE] |>
    as.data.frame() |>
    tibble::rownames_to_column("sample") |>
    pivot_longer(-sample, names_to = "species", values_to = value_name)
}

long_before <- to_long(relabund_before, "value")
long_after  <- to_long(relabund_after,  "value")
long_change <- to_long(change,          "value")

species_levels <- top_species  # already ranked most -> least abundant

fix_levels <- function(df) {
  df$sample  <- factor(df$sample,  levels = sample_levels)
  df$species <- factor(df$species, levels = species_levels)
  df
}
long_before <- fix_levels(long_before)
long_after  <- fix_levels(long_after)
long_change <- fix_levels(long_change)

## ---- Shared color scales ----
## Panels 1 & 2 (abundance) share one sequential scale for direct visual comparison.
abund_max <- max(c(long_before$value, long_after$value))
seq_low  <- "#eef2f6"
seq_high <- "#1c3d5a"

## Panel 3 (change) is diverging, symmetric around zero.
change_lim <- max(abs(long_change$value))
div_neg <- "#b34d3e"   # loss (relative abundance went down after decontam)
div_mid <- "#f3f1ec"
div_pos <- "#2f6f5e"   # gain

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
    axis.text.y = element_text(size = 4.5),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.width = unit(3, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )

p1 <- ggplot(long_before, aes(x = species, y = sample, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = seq_low, high = seq_high, limits = c(0, abund_max),
                       labels = percent_format(accuracy = 1), name = "Rel. abund.") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Before decontamination") +
  base_theme

p2 <- ggplot(long_after, aes(x = species, y = sample, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = seq_low, high = seq_high, limits = c(0, abund_max),
                       labels = percent_format(accuracy = 1), name = "Rel. abund.") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "After decontamination") +
  base_theme

p3 <- ggplot(long_change, aes(x = species, y = sample, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = div_neg, mid = div_mid, high = div_pos, midpoint = 0,
                        limits = c(-change_lim, change_lim),
                        labels = percent_format(accuracy = 1), name = "Change") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Change (after - before)") +
  base_theme

## ---- Combine panels side by side with a shared title/subtitle -----------
## No patchwork/gridExtra/cowplot/ggtext in the scrub:fcbb852 container, so panels are
## laid out by hand with base grid (ggplotGrob + viewport), and each panel's control/blank
## sample labels are recolored red by editing the rendered axis-text grob directly (see
## locate_axis_text_y() below), rather than via ggtext markdown.
title    <- "Effect of SCRuB decontamination on relative abundance"
subtitle <- sprintf("Top %d species by mean relative abundance · %d samples · %d control/blank sample(s) shown in red (no after/change data)",
                     length(top_species), length(sample_levels), length(control_samples))

## Opened before building the grobs: ggplotGrob() needs an active graphics device for
## font-metric calculations, and would otherwise silently open (and leave behind) the
## default Rplots.pdf device.
ragg::agg_png(opt$out, width = 16, height = 9, units = "in", res = 300, background = "white")

grobs <- lapply(list(p1, p2, p3), ggplotGrob)

## Each panel's y-axis tick labels are one grid text grob holding all labels in order
## (inside the "axis-l" cell, a titleGrob wrapping a single text child); locate that text
## grob so its per-label colors can be overwritten, coloring control/blank sample IDs red
## without ggtext.
locate_axis_text_y <- function(axis_grob) {
  axis_cell  <- axis_grob$children$axis
  ti         <- which(vapply(axis_cell$grobs, inherits, logical(1), what = "titleGrob"))[1]
  title_grob <- axis_cell$grobs[[ti]]
  text_name  <- title_grob$childrenOrder[1]
  list(axis_cell = axis_cell, title_idx = ti, title_grob = title_grob,
       text_name = text_name, text_grob = title_grob$children[[text_name]])
}

for (i in seq_along(grobs)) {
  axis_l_idx <- which(grobs[[i]]$layout$name == "axis-l")
  loc <- locate_axis_text_y(grobs[[i]]$grobs[[axis_l_idx]])
  stopifnot(identical(loc$text_grob$label, sample_levels))

  label_colors <- ifelse(is_control[sample_levels], CONTROL_COLOR, loc$text_grob$gp$col)
  loc$text_grob$gp$col <- label_colors
  loc$title_grob$children[[loc$text_name]] <- loc$text_grob
  loc$axis_cell$grobs[[loc$title_idx]] <- loc$title_grob
  grobs[[i]]$grobs[[axis_l_idx]]$children$axis <- loc$axis_cell
}

grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 3, ncol = length(grobs),
  heights = unit(c(1.6, 1.1, 1), c("lines", "lines", "null"))
)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = seq_along(grobs)))
grid.text(title, x = unit(2, "mm"), just = "left", gp = gpar(fontsize = 15, fontface = "bold"))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = seq_along(grobs)))
grid.text(subtitle, x = unit(2, "mm"), just = "left", gp = gpar(fontsize = 10, col = "grey40"))
popViewport()

for (i in seq_along(grobs)) {
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i))
  grid.draw(grobs[[i]])
  popViewport()
}
popViewport()
invisible(dev.off())
