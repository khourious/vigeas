# vigeas-nanopore-depth-coverage
# --> updated: 26 Jun 2023
# ----> Laise de Moraes [https://lpmor22.github.io/]
# ------> Khouri Lab, Gonçalo Moniz Institute, FIOCRUZ, Brazil

if(!interactive()) pdf(NULL)

library("cowplot")
library("dplyr")
library("ggplot2")
library("patchwork")
library("plyr")
library("readr")
library("svglite")

args <- commandArgs(TRUE)
id_sample <- args[3]
primer_scheme <- args[4]
output <- args[5]
print(args)

input_depth <- read.delim(args[1], header = FALSE)
depth_coverage <- data.frame(position = input_depth$V2, depth = input_depth$V3)

ref_seq <- list(
  "ARTIC" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "FIOCRUZ-IOC" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "MIDNIGHT" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ZikaAsian" = "Genome reference: KJ776791.2 (ZIKV)",
  "DENGUESEQ1" = "Genome reference: NC_001477.1 (DENV-1)",
  "DENV1" = "Genome reference: NC_001477.1 (DENV-1)",
  "DENGUESEQ2" = "Genome reference: NC_001474.2 (DENV-2)",
  "DENV2" = "Genome reference: NC_001474.2 (DENV-2)",
  "DENGUESEQ3" = "Genome reference: NC_001475.2 (DENV-3)",
  "DENV3" = "Genome reference: NC_001475.2 (DENV-3)",
  "DENGUESEQ4" = "Genome reference: NC_002640.1 (DENV-4)",
  "DENV4" = "Genome reference: NC_002640.1 (DENV-4)",
  "ChikAsianECSA" = "Genome reference: KP164568.1 (CHIKV)",
  "HTLV1" = "Genome reference: J02029.1 (HTLV-1)",
  "WNV400" = "Genome reference: NC_009942.1 (WNV)")
primer_scheme_2 <- ref_seq[[primer_scheme]]

if (primer_scheme == "ARTIC" || primer_scheme == "FIOCRUZ-IOC" || primer_scheme == "MIDNIGHT") {
  # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
  # https://doi.org/10.1038/s41586-020-2286-9
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 265,
                  "ORFs", "ORF1a", 266, 13468,
                  "ORFs", "ORF1b", 13468, 21555,
                  "Structural proteins", "S", 21563, 25384,
                  "Structural proteins", "E", 26245, 26472,
                  "Structural proteins", "M", 26523, 27191,
                  "Structural proteins", "N", 28274, 29533,
                  "UTR", "3'UTR", 29675, 29903)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "NSPs", "nsp1", 266, 805,
                  "NSPs", "nsp2", 806, 2719,
                  "NSPs", "nsp3", 2720, 8554,
                  "NSPs", "nsp4", 8555, 10054,
                  "NSPs", "nsp5", 10055, 10972,
                  "NSPs", "nsp6", 10973, 11842,
                  "NSPs", "nsp8", 12092, 12685,
                  "NSPs", "nsp10", 13025, 13441,
                  "NSPs", "nsp12", 13442, 16236,
                  "NSPs", "nsp13", 16237, 18039,
                  "NSPs", "nsp14", 18040, 19620,
                  "NSPs", "nsp15", 19621, 20658,
                  "NSPs", "nsp16", 20659, 21552,
                  "Accessory factors", "3ab", 25393, 26220,
                  "Accessory factors", "6", 27202, 27387,
                  "Accessory factors", "7a", 27394, 27759,
                  "Accessory factors", "7b", 27756, 27887,
                  "Accessory factors", "8", 27894, 28259,
                  "Accessory factors", "9b", 28284, 28577,
                  "Accessory factors", "10", 29558, 28674)
  map3 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "NSPs", "nsp7", 11843, 12091,
                  "NSPs", "nsp9", 12686, 13024,
                  "NSPs", "nsp11", 13442, 13480)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#2F67CD",
      "ORFs" = "#FE0B12",
      "Structural proteins" = "#11961B"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "NSPs" = "#FF8A8E",
      "Accessory factors" = "#D860CF"))
  map3plot <- map3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "NSPs" = "#FF8A8E"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30000)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".sars2-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot / map3plot + plot_layout(nrow = 4, heights = c(3, .4, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                         expand = expansion(0, 0), limits = c(0, 30000)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".sars2-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot / map3plot + plot_layout(nrow = 4, heights = c(3, .4, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "ZikaAsian") {
  # https://www.ncbi.nlm.nih.gov/nuccore/KJ776791.2
  # https://doi.org/10.1371/journal.ppat.1006528
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 107,
                  "Structural proteins", "pr", 474, 752,
                  "Structural proteins", "M", 753, 977,
                  "Structural proteins", "E", 978, 2489,
                  "Non-structural proteins", "NS1", 2490, 3545,
                  "Non-structural proteins", "NS2A", 3546, 4223,
                  "Non-structural proteins", "NS3", 4614, 6464,
                  "Non-structural proteins", "NS4A", 6465, 6845,
                  "Non-structural proteins", "NS4B", 6915, 7667,
                  "Non-structural proteins", "NS5", 7668, 10379,
                  "UTR", "3'UTR", 10380, 10807)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 108, 473,
                  "Non-structural proteins", "NS2B", 4224, 4613,
                  "ORF", "2K", 6846, 6914)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10810)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#FE0B12",
      "Structural proteins" = "#11961B",
      "Non-structural proteins" = "#2F67CD"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10810)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#11961B",
      "Non-structural proteins" = "#2F67CD",
      "ORF" = "#FE0B12"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10807),
                       expand = expansion(0, 0), limits = c(0, 10810)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".zikv-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10807),
                         expand = expansion(0, 0), limits = c(0, 10810)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".zikv-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "DENGUESEQ1" || primer_scheme == "DENV1") {
  # https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 94,
                  "Structural proteins", "prM", 437, 934,
                  "Structural proteins", "E", 935, 2419,
                  "Non-structural proteins", "NS1", 2420, 3475,
                  "Non-structural proteins", "NS2A", 3476, 4129,
                  "Non-structural proteins", "NS3", 4520, 6376,
                  "Non-structural proteins", "NS4A", 6377, 6757,
                  "Non-structural proteins", "NS4B", 6827, 7573,
                  "Non-structural proteins", "NS5", 7574, 10270,
                  "UTR", "3'UTR", 10274, 10735)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 95, 436,
                  "Non-structural proteins", "NS2B", 4130, 4519,
                  "ORF", "2K", 6758, 6826)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#0247fe",
      "Structural proteins" = "#80a3fe",
      "Non-structural proteins" = "#ccdaff"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#80a3fe",
      "Non-structural proteins" = "#ccdaff",
      "ORF" = "#0247fe"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10735),
                       expand = expansion(0, 0), limits = c(0, 10810)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".denv1-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10735),
                         expand = expansion(0, 0), limits = c(0, 10750)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".denv1-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
}}

if (primer_scheme == "DENGUESEQ2" || primer_scheme == "DENV2") {
  # https://www.ncbi.nlm.nih.gov/nuccore/NC_001474.2
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 96,
                  "Structural proteins", "prM", 439, 936,
                  "Structural proteins", "E", 937, 2421,
                  "Non-structural proteins", "NS1", 2422, 3477,
                  "Non-structural proteins", "NS2A", 3478, 4131,
                  "Non-structural proteins", "NS3", 4522, 6375,
                  "Non-structural proteins", "NS4A", 6376, 6756,
                  "Non-structural proteins", "NS4B", 6826, 7569,
                  "Non-structural proteins", "NS5", 7570, 10269,
                  "UTR", "3'UTR", 10273, 10723)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 97, 438,
                  "Non-structural proteins", "NS2B", 4132, 4521,
                  "ORF", "2K", 6757, 6825)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#f23513",
      "Structural proteins" = "#ff9c80",
      "Non-structural proteins" = "#ffded3"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#ff9c80",
      "Non-structural proteins" = "#ffded3",
      "ORF" = "#f23513"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".denv2-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10735),
                         expand = expansion(0, 0), limits = c(0, 10750)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".denv2-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "DENGUESEQ3" || primer_scheme == "DENV3") {
  # https://www.ncbi.nlm.nih.gov/nuccore/NC_001475.2
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 94,
                  "Structural proteins", "prM", 437, 934,
                  "Structural proteins", "E", 935, 2413,
                  "Non-structural proteins", "NS1", 2414, 3469,
                  "Non-structural proteins", "NS2A", 3470, 4123,
                  "Non-structural proteins", "NS3", 4514, 6370,
                  "Non-structural proteins", "NS4A", 6371, 6751,
                  "Non-structural proteins", "NS4B", 6821, 7564,
                  "Non-structural proteins", "NS5", 7565, 10264,
                  "UTR", "3'UTR", 10268, 10707)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 95, 436,
                  "Non-structural proteins", "NS2B", 4124, 4513,
                  "ORF", "2K", 6752, 6820)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#c101ba",
      "Structural proteins" = "#fe57f8",
      "Non-structural proteins" = "#ffc7fd"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#fe57f8",
      "Non-structural proteins" = "#ffc7fd",
      "ORF" = "#c101ba"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".denv3-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10707),
                         expand = expansion(0, 0), limits = c(0, 10750)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".denv3-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "DENGUESEQ4" || primer_scheme == "DENV4") {
  # https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 101,
                  "Structural proteins", "prM", 441, 938,
                  "Structural proteins", "E", 939, 2423,
                  "Non-structural proteins", "NS1", 2424, 3479,
                  "Non-structural proteins", "NS2A", 3480, 4133,
                  "Non-structural proteins", "NS3", 4524, 6377,
                  "Non-structural proteins", "NS4A", 6378, 6758,
                  "Non-structural proteins", "NS4B", 6828, 7562,
                  "Non-structural proteins", "NS5", 7563, 10262,
                  "UTR", "3'UTR", 10266, 10649)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 102, 440,
                  "Non-structural proteins", "NS2B", 4134, 4523,
                  "ORF", "2K", 6759, 6827)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#158e01",
      "Structural proteins" = "#24f901",
      "Non-structural proteins" = "#b1ffa4"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#24f901",
      "Non-structural proteins" = "#b1ffa4",
      "ORF" = "#158e01"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10649),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".denv4-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10707),
                         expand = expansion(0, 0), limits = c(0, 10750)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".denv4-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "ChikAsianECSA") {
  # https://www.ncbi.nlm.nih.gov/nuccore/KP164568.1
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 77,
                  "Non-structural proteins", "NS1", 78, 1682,
                  "Non-structural proteins", "NS2", 1683, 4076,
                  "Non-structural proteins", "NS3", 4077, 5648,
                  "Non-structural proteins", "NS4", 5667, 7499,
                  "UTR", "3'UTR", 11315, 11812)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 7568, 8359,
                  "Structural proteins", "E3", 8360, 8542,
                  "Structural proteins", "E2", 8543, 9811,
                  "Structural proteins", "6K", 9812, 9994,
                  "Structural proteins", "E1", 9995, 11311)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 12000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#FE0B12",
      "Structural proteins" = "#11961B",
      "Non-structural proteins" = "#2F67CD"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 12000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#11961B",
      "Non-structural proteins" = "#2F67CD",
      "ORF" = "#FE0B12"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 11812),
                       expand = expansion(0, 0), limits = c(0, 12000)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".chikv-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 11812),
                         expand = expansion(0, 0), limits = c(0, 12000)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".chikv-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "HTLV1") {
  # https://www.ncbi.nlm.nih.gov/nuccore/J02029.1
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "gag", "gag", 824, 2113,
                  "pol", "pol", 2114, 5210,
                  "pX", "pX", 6857, 8382)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "LTR", "5'LTR", 23, 777,
                  "pro", "pro", 1960, 2778,
                  "env", "env", 5203, 6669,
                  "LTR", "3'LTR", 8301, 9055)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "LTR" = "#5A5A5A",
      "gag" = "#FE0B12",
      "pol" = "#2F67CD",
      "pX" = "#FE810F"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "env" = "#11961B",
      "pro" = "#2F67CD"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 2000, 4000, 6000, 8000, 9068),
                       expand = expansion(0, 0), limits = c(0, 9100)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".htlv1-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 2000, 4000, 6000, 8000, 9068),
                         expand = expansion(0, 0), limits = c(0, 12000)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".htlv1-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}

if (primer_scheme == "WNV400") {
  # https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1
  map1 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "UTR", "5'UTR", 1, 96,
                  "Structural proteins", "pr", 466, 741,
                  "Structural proteins", "M", 742, 966,
                  "Structural proteins", "E", 967, 2469,
                  "Non-structural proteins", "NS1", 2470, 3525,
                  "Non-structural proteins", "NS2A", 3526, 4218,
                  "Non-structural proteins", "NS3", 4612, 6468,
                  "Non-structural proteins", "NS4A", 6469, 6846,
                  "Non-structural proteins", "NS4B", 6916, 7680,
                  "Non-structural proteins", "NS5", 7681, 10395,
                  "UTR", "3'UTR", 10399, 11029)
  map2 <- tribble(~"class", ~"gene", ~"start", ~"end",
                  "Structural proteins", "C", 97, 411,
                  "Non-structural proteins", "NS2B", 4219, 4611,
                  "ORF", "2K", 6847, 6915)
  map1plot <- map1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 11060)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "UTR" = "#FE0B12",
      "Structural proteins" = "#11961B",
      "Non-structural proteins" = "#2F67CD"))
  map2plot <- map2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
              linewidth = .2, colour = "#000000", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 11060)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c(
      "Structural proteins" = "#11961B",
      "Non-structural proteins" = "#2F67CD",
      "ORF" = "#FE0B12"))
  depcov1 <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11029),
                       expand = expansion(0, 0), limits = c(0, 11060)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme_light(base_size = 10) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(hjust = 1, size = 8)) +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
  output1 <-  paste0(output, ".wnv-coverage.pdf")
  plot1 <- depcov1 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
  save_plot(output1, plot1, base_height = 5, base_width = 16)
  if (file.size(args[2]) > 0) {
    contamination_bed <- read.delim(args[2], header = FALSE)
    contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)
    depcov2 <- ggplot() +
      geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
      geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
      labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11029),
                         expand = expansion(0, 0), limits = c(0, 11060)) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(hjust = 1, size = 8)) +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#5A5A5A")
    output2 <-  paste0(output, ".wnv-coverage.contamination.pdf")
    plot2 <- depcov2 / map1plot / map2plot + plot_layout(nrow = 3, heights = c(3, .3, .3))
    save_plot(output2, plot2, base_height = 5, base_width = 16)
  }}