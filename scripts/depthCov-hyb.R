if (!require("pacman")) {install.packages("pacman", dependencies = TRUE)}
pacman::p_load("cowplot",
               "dplyr",
               "ggplot2",
               "patchwork",
               "plyr",
               "qpdf",
               "readr",
               "tidyverse",
               "rstudioapi",
               "svglite")

path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

input <- list.files(pattern = "\\.depth\\.tsv$",
                    full.names = TRUE,
                    recursive = TRUE)

ref_seq <- list(
  "K03455" =
    "Genome reference: K03455.1 (Lentivirus humimdef1 - HIV-1)",
  "U63715" =
    "Genome reference: U63715.1 (Pegivirus hominis 1 - HPgV-1)",
  "AF345523" =
    "Genome reference: AF345523.1 (Alphatorquevirus homin5 - TTV5)",
  "AB025946" =
    "Genome reference: AB025946.2 (Alphatorquevirus homin19 - TTV19)",
  "AB038621" =
    "Genome reference: AB038621.1 (Alphatorquevirus homin29 - TTV29)",
  "JN555585" =
    "Genome reference: JN555585.1 (Simplexvirus humanalpha1 - HSV-1)",
  "JN561323" =
    "Genome reference: JN561323.2 (Simplexvirus humanalpha2 - HSV-2)",
  "AF484424" =
    "Genome reference: AF484424.1 (Orthobunyavirus oropoucheense, L segment - OROV)",
  "AF441119" =
    "Genome reference: AF441119.1 (Orthobunyavirus oropoucheense, M segment - OROV)",
  "AY237111" =
    "Genome reference: AY237111.1 (Orthobunyavirus oropoucheense, S segment - OROV)")

for (i in input) {
  id_sample <- strsplit(basename(i), "\\.")[[1]][1]
  target <- strsplit(basename(i), "\\.")[[1]][2]
  target_refseq <- ref_seq[[target]]
  
  if (file.info(i)$size == 0) {next}
  
  df <- read.delim(i, sep = "\t", header = FALSE)
  if (nrow(df) == 0) {next}
  
  depth_coverage <- data.frame(position = df$V2, depth = df$V3)
  
  if (target == "K03455") {
    # https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .4, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 634, 790, 2292, 5041, 5619, 5831,
                                    6045, 6062, 6310, 8379, 8424, 8653,
                                    8797, 9086, 9417, 9719, 9728),
                         expand = expansion(0, 0), limits = c(0, 9800),
                         guide = guide_axis(n.dodge = 2)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    map2genome1 <- tribble(~"gene", ~"start", ~"end",
                           "5'LTR", 1, 634,
                           "gag", 790, 2292,
                           "vif", 5041, 5619,
                           "tat", 8379, 8424,
                           "nef", 8797, 9417)
    map2plot1 <- map2genome1 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "red", colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")
    map2genome2 <- tribble(~"gene", ~"start", ~"end",
                           "tat", 5831, 6045,
                           "vpu", 6062, 6310,
                           "rev", 8379, 8653,
                           "3'LTR", 9086, 9719)
    map2plot2 <- map2genome2 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "red", colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")
    map2genome3 <- tribble(~"gene", ~"start", ~"end",
                           "pol", 2085, 5096,
                           "vpr", 5559, 5850,
                           "rev", 5970, 6045,
                           "env", 6225, 8795)
    map2plot3 <- map2genome3 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "red", colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")  
    id_sample <- paste0(id_sample, ".HIV1-coverage.pdf")
    plot <- depcov / map2plot1 / map2plot2 / map2plot3 +
      plot_layout(nrow = 4, heights = c(3, .3, .3, .3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
  
  if (target == "U63715") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .8, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 2000, 3000, 4000, 5000,
                                    6000, 7000, 8000, 9000, 9367),
                         expand = expansion(0, 0), limits = c(0, 9400)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".HPgV1-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
  
  if (target == "AF345523") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .8, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 500, 1000, 1500, 2000, 2500,
                                    3000, 3229),
                         expand = expansion(0, 0), limits = c(0, 3250),
                         guide = guide_axis(n.dodge = 2)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".TTV5-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
  
  if (target == "AB025946") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .8, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 500, 1000, 1500, 2000, 2500,
                                    3000, 3500, 3808),
                         expand = expansion(0, 0), limits = c(0, 3850),
                         guide = guide_axis(n.dodge = 2)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".TTV19-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
  
  if (target == "AB038621") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .8, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 500, 1000, 1500, 2000, 2500,
                                    3000, 3500, 3676),
                         expand = expansion(0, 0), limits = c(0, 3700),
                         guide = guide_axis(n.dodge = 2)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".TTV29-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
  
  if (target == "JN555585") {
    # https://doi.org/10.3389/fimmu.2018.02099
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .4, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 9213, 9338, 116926, 118777, 127151,
                                    127173, 131460, 132101, 146096, 146737,
                                    151023, 152165, 152222),
                         expand = expansion(0, 0), limits = c(0, 154000),
                         guide = guide_axis(n.dodge = 2)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    map2genome1 <- tribble(~"gene", ~"start", ~"end",
                           "TR[L]", 1, 9213,
                           "TR[S]", 146737, 151023)
    map2plot1 <- map2genome1 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "green",
                colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene),
                parse = TRUE, size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 154000)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")
    map2genome2 <- tribble(~"gene", ~"start", ~"end",
                           "UL[S]", 9338, 116926,
                           "U[S]", 132101, 146096)
    map2plot2 <- map2genome2 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "yellow",
                colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene),
                parse = TRUE, size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 154000)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")
    map2genome3 <- tribble(~"gene", ~"start", ~"end",
                           "IR[L]", 118777, 127151,
                           "IR[S]", 127173, 131460)
    map2plot3 <- map2genome3 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "red",
                colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene),
                parse = TRUE, size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 154000)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")  
    id_sample <- paste0(id_sample, ".HSV1-coverage.pdf")
    plot <- depcov / map2plot1 / map2plot2 / map2plot3 +
      plot_layout(nrow = 4, heights = c(3, .3, .3, .3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
  
  if (target == "JN561323") {
    # https://doi.org/10.1128/jvi.72.3.2010-2021.1998
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .4, colour = "black") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 9297, 9298, 117782, 117987, 127249,
                                    128079, 132898, 133707, 148035, 148036,
                                    154746, 154675),
                         expand = expansion(0, 0), limits = c(0, 156000),
                         guide = guide_axis(n.dodge = 2)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    map2genome1 <- tribble(~"gene", ~"start", ~"end",
                           "TR[L]", 1, 9297,
                           "TR[S]", 148036, 154746)
    map2plot1 <- map2genome1 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "green",
                colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene),
                parse = TRUE, size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 156000)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")
    map2genome2 <- tribble(~"gene", ~"start", ~"end",
                           "UL[S]", 9298, 117782,
                           "U[S]", 133707, 148035)
    map2plot2 <- map2genome2 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "yellow",
                colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene),
                parse = TRUE, size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 156000)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")
    map2genome3 <- tribble(~"gene", ~"start", ~"end",
                           "IR[L]", 117987, 127249,
                           "IR[S]", 128079, 132898)
    map2plot3 <- map2genome3 %>% ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
                linewidth = .2, fill = "red",
                colour = "darkgray", alpha = .3) +
      geom_text(aes(x = (start + end) / 2, y = 9, label = gene),
                parse = TRUE, size = 6) +
      scale_x_continuous(expand = expansion(0, 0), limits = c(0, 156000)) +
      theme_void() + theme(legend.position = "none") +
      coord_cartesian(clip = "off")  
    id_sample <- paste0(id_sample, ".HSV2-coverage.pdf")
    plot <- depcov / map2plot1 / map2plot2 / map2plot3 +
      plot_layout(nrow = 4, heights = c(3, .3, .3, .3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }

  if (target == "AF484424") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .6, colour = "darkblue") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 1000, 2000, 3000, 4000, 5000,
                                    6000, 6846),
                         expand = expansion(0, 0), limits = c(0, 6900)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".OROVL-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }

  if (target == "AF441119") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .6, colour = "darkred") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 500, 1000, 1500, 2000, 2500, 3000,
                                    3500, 4000, 4500, 4385),
                         expand = expansion(0, 0), limits = c(0, 4400)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".OROVM-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }

  if (target == "AY237111") {
    depcov <- ggplot() +
      geom_line(data = depth_coverage, aes(x = position, y = depth),
                linewidth = .6, colour = "darkgreen") +
      labs(title = paste0(id_sample), subtitle = paste0(target_refseq),
           y = "Per base coverage (x)", x = NULL) +
      scale_x_continuous(breaks = c(1, 100, 200, 300, 400, 500, 600, 700, 754),
                         expand = expansion(0, 0), limits = c(0, 800)) +
      theme_light(base_size = 10) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                    expand = expansion(0, 0)) +
      theme(panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title.y = element_text(angle = 90, size = 14),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(hjust = 1, size = 9)) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "black")
    id_sample <- paste0(id_sample, ".OROVS-coverage.pdf")
    plot <- depcov + plot_layout(nrow = 1, heights = c(3))
    save_plot(id_sample, plot, base_height = 7, base_width = 20)
  }
}

pdf <- list.files(pattern = "\\coverage.pdf$", full.names = TRUE)

if(length(pdf) > 0) {
  combined <- "CoverageDepth.pdf"
  pdf_combine(pdf, combined)
}
