#!/usr/bin/env Rscript
#
# 12_figures.R
# ------------
# The four figures for the paper and the README.
#
#   Fig 1  Sampling: Bangladeshi sequences per year, by lineage and sequence type.
#          The honest picture of what data exists, including the gaps.
#   Fig 2  Whole-genome ML tree: tips by origin, support on the branches that
#          carry the argument.
#   Fig 3  Dated tree: the lineage running from ~2016 to 2025, on a calendar axis.
#   Fig 4  Introductions: the same dated tree with branches coloured by the
#          inferred region, so the two introductions are visible.
#
# Each figure is wrapped in tryCatch and runs independently, so one failure does
# not lose the others. The script reports which succeeded at the end.
#
# Outputs: results/figures/fig1_sampling.{png,pdf} ... fig4_introductions.{png,pdf}
#          results/figures/sessionInfo.txt   (for the methods section)
#
# Run from the project folder:
#   Rscript scripts/12_figures.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggtree)
  library(treeio)
  library(patchwork)
})

# Which results directory to read. Defaults to results/ (the step-by-step run);
# pass results_pipeline/ to draw the figures from the Nextflow output instead.
#   Rscript scripts/12_figures.R                     # results/
#   Rscript scripts/12_figures.R results_pipeline    # the pipeline's output
.args <- commandArgs(trailingOnly = TRUE)
RES <- if (length(.args) >= 1 && nzchar(.args[1])) .args[1] else "results"
if (!dir.exists(RES)) stop("results directory not found: ", RES)
message("reading results from: ", RES)

# TreeTime's time-scaled trees contain a few very small NEGATIVE branch lengths,
# a rounding artefact of converting divergence to time. ggtree's default is to
# drop those edges, which silently distorts the tree. Keeping them is correct:
# they are effectively zero-length, not real.
options(ignore.negative.edge = TRUE)

OUT <- file.path(RES, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# small helper so every path below is relative to RES
R <- function(...) file.path(RES, ...)

# A single palette used across every figure, so colours mean the same thing
# everywhere. Bangladesh is the only saturated colour; everything else is muted,
# which is what makes the figures readable at a glance.
PAL <- c(
  "Bangladesh"         = "#C1272D",
  "South Asia (other)" = "#E09F3E",
  "Southeast Asia"     = "#335C67",
  "East Asia"          = "#9E7BB5",
  "elsewhere"          = "#9A9A9A",
  "unknown"            = "#D6D6D6"
)
LINEAGE_PAL <- c(
  "3I_A.2"  = "#C1272D",
  "3III_A"  = "#335C67",
  "3II_B"   = "#E09F3E",
  "other"   = "#9A9A9A"
)

theme_pub <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(colour = "grey25", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey25", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = base + 2),
      plot.subtitle = element_text(colour = "grey35", size = base - 1),
      plot.caption = element_text(colour = "grey45", size = base - 2, hjust = 0),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# Long subtitles and captions were being cut off at the plot edge.
wrap <- function(x, width = 95) paste(strwrap(x, width = width), collapse = "\n")

# Tips TreeTime's clock filter rejected. Their dates in the time-scaled tree are
# fabricated (2075-2086 for 2024 samples), and because ggtree anchors the time
# axis to the most distant tip, leaving them in shifted every real tip ~25 years
# into the past. They are excluded from rate estimation anyway, so they are
# dropped from the dated figures and the count is reported in the caption.
clock_outlier_tips <- function() {
  f <- R("dating", "wg", "auspice_tree.json")
  if (!file.exists(f)) {
    message("  no auspice_tree.json at ", f, " - cannot identify clock outliers")
    return(character(0))
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("  jsonlite not installed - cannot identify clock outliers. ",
            "install.packages('jsonlite')")
    return(character(0))
  }
  j <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  out <- character(0)
  walk <- function(n) {
    if (is.null(n$children) || !length(n$children)) {
      bb <- n$node_attrs$bad_branch$value
      if (!is.null(bb) && identical(bb, "Yes")) out <<- c(out, n$name)
    } else for (k in n$children) walk(k)
  }
  walk(j$tree)
  unique(out)
}

drop_outliers <- function(tr) {
  bad <- clock_outlier_tips()
  labs <- if (inherits(tr, "treedata")) tr@phylo$tip.label else tr$tip.label
  hit <- intersect(bad, labs)
  if (!length(hit)) return(list(tree = tr, n = 0))
  message("  dropping ", length(hit), " clock-filtered tip(s) with fabricated ",
          "dates: ", paste(hit, collapse = ", "))
  # ape::drop.tip dispatches on 'phylo' only; treeio supplies the method for
  # 'treedata' and keeps the node annotations Fig 4 needs.
  res <- try(
    if (inherits(tr, "treedata")) treeio::drop.tip(tr, hit) else ape::drop.tip(tr, hit),
    silent = TRUE)
  if (inherits(res, "try-error")) {
    message("    could not drop them (", trimws(as.character(res)), ")")
    message("    continuing with all tips - the axis below is still correct, ",
            "because it is anchored on the root date, not on the furthest tip")
    return(list(tree = tr, n = 0))
  }
  list(tree = res, n = length(hit))
}

# The date of the root, read from the auspice tree. Branch lengths in TreeTime's
# time tree are in YEARS, so calendar date = root date + distance from root.
# This is exact, and unlike ggtree's mrsd it does not depend on which tip is
# furthest from the root - which is what made the first version wrong by 25
# years, because clock-filtered tips carry fabricated dates out to 2086.
root_num_date <- function() {
  f <- R("dating", "wg", "auspice_tree.json")
  if (!file.exists(f) || !requireNamespace("jsonlite", quietly = TRUE)) return(NA_real_)
  j <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  v <- j$tree$node_attrs$num_date$value
  if (is.null(v)) NA_real_ else as.numeric(v)
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  message("  wrote ", name, ".png and .pdf")
}

ok <- character(); failed <- character()
attempt <- function(label, fn) {
  message("\n== ", label, " ==")
  tryCatch({ fn(); ok <<- c(ok, label) },
           error = function(e) {
             message("  FAILED: ", conditionMessage(e))
             failed <<- c(failed, label)
           })
}

# ---------------------------------------------------------------- helpers ----

read_meta <- function() {
  read_tsv(R("metadata", "metadata_clades.tsv"), show_col_types = FALSE)
}

# Node labels from IQ-TREE look like "99.9/100" = SH-aLRT / UFBoot.
split_support <- function(lab) {
  lab <- as.character(lab)
  has <- grepl("/", lab, fixed = TRUE)
  alrt <- rep(NA_real_, length(lab)); boot <- rep(NA_real_, length(lab))
  parts <- strsplit(lab[has], "/", fixed = TRUE)
  alrt[has] <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 1)))
  boot[has] <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
  list(alrt = alrt, boot = boot)
}

read_tree_safe <- function(path) {
  # treeio::read.beast handles [&key=value] annotations; fall back to plain nexus
  t <- try(treeio::read.beast(path), silent = TRUE)
  if (inherits(t, "try-error")) t <- treeio::read.nexus(path)
  t
}

tip_table <- function(stem) {
  traits <- read_tsv(R("datasets", sprintf("%s_traits.tsv", stem)),
                     show_col_types = FALSE)
  dates <- read_tsv(R("datasets", sprintf("%s_dates.tsv", stem)),
                    show_col_types = FALSE)
  traits %>%
    left_join(dates, by = "name") %>%
    mutate(
      year = suppressWarnings(as.integer(substr(as.character(date), 1, 4))),
      origin = ifelse(region %in% names(PAL), region, "elsewhere"),
      origin = factor(origin, levels = names(PAL))
    )
}

# ------------------------------------------------- Fig 1: what data exists ----

attempt("Fig 1 - sampling", function() {
  md <- read_meta()
  bd <- md %>%
    filter(country == "Bangladesh", !is.na(year)) %>%
    mutate(
      lineage = ifelse(clade %in% names(LINEAGE_PAL), clade, "other"),
      lineage = factor(lineage, levels = names(LINEAGE_PAL)),
      kind = case_when(
        length >= 10000 ~ "whole genome",
        length >= 1400  ~ "E-gene length",
        TRUE            ~ "shorter fragment"
      ),
      kind = factor(kind, levels = c("whole genome", "E-gene length",
                                     "shorter fragment"))
    )

  yrs <- seq(min(bd$year), max(bd$year))
  missing_yrs <- setdiff(yrs, unique(bd$year))
  # %in% rather than ==: md$country contains NAs, and `== "Bangladesh"`
  # yields NA for those, which propagates through sum() and printed "NA"
  n_undated <- sum(md$country %in% "Bangladesh" & is.na(md$year))

  p1 <- ggplot(bd, aes(x = year, fill = lineage)) +
    geom_bar(width = 0.78) +
    scale_fill_manual(values = LINEAGE_PAL, drop = FALSE) +
    scale_x_continuous(breaks = yrs) +
    labs(title = "DENV-3 sequences from Bangladesh, by year and lineage",
         subtitle = wrap(sprintf(
           "%d dated sequences (a further %d carry no collection date). Years with no public data: %s",
           nrow(bd), n_undated,
           if (length(missing_yrs)) paste(missing_yrs, collapse = ", ") else "none")),
         x = NULL, y = "sequences") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(bd, aes(x = year, fill = kind)) +
    geom_bar(width = 0.78) +
    scale_fill_manual(values = c("whole genome" = "#335C67",
                                 "E-gene length" = "#7FA8B0",
                                 "shorter fragment" = "#D6D6D6")) +
    scale_x_continuous(breaks = yrs) +
    labs(subtitle = "The same sequences by how much of the genome they cover",
         x = NULL, y = "sequences",
         caption = paste(
           wrap(paste("Whole genomes are what the dated analysis uses;",
                      "E-gene-length sequences carry the year-to-year continuity.")),
           "\nSource: NCBI Virus, taxid 11069, downloaded 31 October 2026.")) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_fig(p1 / p2, "fig1_sampling", 9, 8)
})

# ------------------------------------------ Fig 2: whole-genome ML tree -------

attempt("Fig 2 - whole-genome tree", function() {
  # IQ-TREE's .treefile is plain Newick, NOT nexus. read.beast cannot parse it
  # and the nexus fallback leaves node.label empty, which surfaced as
  # "argument of length 0" and a max(ii) warning. ape::read.tree reads it
  # correctly and keeps the "SH-aLRT/UFBoot" labels as strings.
  tr <- ape::read.tree(R("trees", "wg.treefile"))
  if (is.null(tr$node.label) || !length(tr$node.label)) {
    stop("no internal node labels in wg.treefile - IQ-TREE should write ",
         "'SH-aLRT/UFBoot' labels when run with -B and --alrt")
  }
  tips <- tip_table("wg")

  # Keep only well-supported labels, written straight into node.label. Doing it
  # this way avoids a second %<+% join, which is the fragile part of ggtree.
  sup <- split_support(tr$node.label)
  solid <- !is.na(sup$alrt) & !is.na(sup$boot) & sup$alrt >= 80 & sup$boot >= 95
  tr$node.label <- ifelse(solid,
                          sprintf("%.0f/%.0f", sup$alrt, sup$boot), "")
  message("  well-supported internal branches: ", sum(solid), " of ", length(solid))

  p <- ggtree(tr, linewidth = 0.35, colour = "grey30") %<+%
    (tips %>% select(label = name, origin, clade, year)) +
    geom_tippoint(aes(colour = origin), size = 1.5) +
    scale_colour_manual(values = PAL, drop = FALSE, na.value = "grey80") +
    # geom_treescale takes 'color', not 'colour', and the spelling is not
    # interchangeable here as it is in ggplot2
    geom_treescale(width = 0.01, fontsize = 3, color = "grey40") +
    labs(title = "Whole-genome maximum-likelihood tree",
         subtitle = wrap(paste(
           "DENV-3 genotype I, plus the 2021 Bangladeshi genotype III (3III_A)",
           "detection and a 3II_B outgroup. Tip colour = where the sequence was",
           "collected. Labels show SH-aLRT / ultrafast bootstrap, on",
           "well-supported branches only (>=80 and >=95).")),
         caption = paste(
           wrap(paste("IQ-TREE 3.1.4, model GTR+F+I+G4, 1000 ultrafast bootstrap",
                      "and 1000 SH-aLRT replicates.")),
           wrap(paste("Rooted by molecular clock: lineage 3I is a paraphyletic",
                      "parent label and cannot serve as an outgroup.")),
           sep = "\n")) +
    theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), panel.grid = element_blank())

  # check_overlap stops the dense parts of the tree turning into a smear of
  # overlapping support labels
  p <- p + geom_nodelab(size = 2.0, hjust = 1.15, colour = "grey35",
                        check_overlap = TRUE)
  save_fig(p, "fig2_wg_tree", 9, 11)
})

# ------------------------------------------------- Fig 3: dated tree ---------

attempt("Fig 3 - dated tree", function() {
  tr <- read_tree_safe(R("dating", "wg", "timetree.nexus"))
  dropped <- drop_outliers(tr); tr <- dropped$tree
  tips <- tip_table("wg")
  root <- root_num_date()
  if (is.na(root)) stop("could not read the root date from auspice_tree.json")
  message("  root date: ", round(root, 2), " - axis is root + branch length")

  p <- ggtree(tr, linewidth = 0.35, colour = "grey30") %<+%
    (tips %>% select(label = name, origin, clade, year)) +
    theme_tree2() +
    scale_x_continuous(labels = function(b) round(b + root)) +
    geom_tippoint(aes(colour = origin), size = 1.5) +
    scale_colour_manual(values = PAL, drop = FALSE, na.value = "grey80") +
    labs(title = "Time-scaled tree: one lineage, continuously present",
         subtitle = wrap(paste("Most recent common ancestor of the Bangladeshi",
                               "3I_A.2 sequences: 2016.3 (95% CI 2015.8-2016.8)")),
         caption = paste(
           wrap(paste("TreeTime 0.12.1, clock-rooted, --covariation. Rate",
                      "6.81e-4 substitutions/site/year (+/-1.97e-5), R^2 0.94.")),
           wrap(sprintf(paste("PQ815901.1 excluded as a probable recombinant or",
                              "chimeric assembly. %d clock-filtered tip(s) omitted:",
                              "their inferred dates are not real. No Bangladeshi",
                              "sequence is a clock outlier."), dropped$n)),
           sep = "\n")) +
    theme_pub() +
    theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  save_fig(p, "fig3_dated_tree", 9, 11)
})

# --------------------------------------------- Fig 4: introductions ----------

attempt("Fig 4 - introductions", function() {
  path <- R("dating", "wg", "mugration_region", "annotated_tree.nexus")
  tr <- read_tree_safe(path)
  dropped <- drop_outliers(tr); tr <- dropped$tree
  tips <- tip_table("wg")
  root <- root_num_date()
  if (is.na(root)) stop("could not read the root date from auspice_tree.json")

  dat <- as_tibble(tr)
  region_col <- intersect(c("region", "Region"), names(dat))
  if (!length(region_col)) {
    stop(paste("no 'region' annotation found in", path,
               "- columns present:", paste(names(dat), collapse = ", ")))
  }
  dat$region_clean <- ifelse(dat[[region_col[1]]] %in% names(PAL),
                             dat[[region_col[1]]], "elsewhere")

  p <- (ggtree(tr, linewidth = 0.5) %<+%
          (dat %>% select(node, region_clean))) +
    aes(colour = region_clean) +
    theme_tree2() +
    scale_x_continuous(labels = function(b) round(b + root)) +
    scale_colour_manual(values = PAL, drop = FALSE, na.value = "grey80") +
    geom_tippoint(size = 1.3) +
    labs(title = "Inferred location along the tree: two introductions",
         subtitle = wrap(paste(
           "Branch colour = region inferred by TreeTime mugration. One",
           "introduction carries every Bangladeshi 3I_A.2 sequence; a second is",
           "the 2021 3III_A detection.")),
         caption = paste(
           wrap(paste("The NUMBER of introductions is consistent between the",
                      "whole-genome and E-gene datasets (two in both).")),
           wrap(paste("The geographic SOURCE is not: region-level inference gives",
                      "East Asia from whole genomes and Southeast Asia from the E",
                      "gene, so no source is claimed.")),
           sep = "\n")) +
    theme_pub() +
    theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  save_fig(p, "fig4_introductions", 9, 11)
})

# ---------------------------------------------------------------- wrap up ----

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

message("\n================ SUMMARY ================")
for (x in ok) message("  OK      ", x)
for (x in failed) message("  FAILED  ", x)
message("\nFigures in ", OUT, "/   (R package versions in sessionInfo.txt)")
if (length(failed)) {
  message("\nFor any failure above, paste the message - it names the column or",
          "\nfile that did not match, which is what needs fixing.")
}
