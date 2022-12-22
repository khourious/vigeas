# vigeas-depth-coverage
# --> updated: 21 Dec 2022
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
output <- args[4]
print(args)

input_depth <- read.delim(args[1], header = FALSE)
depth_coverage <- data.frame(position = input_depth$V2, depth = input_depth$V3)

contamination_bed <- read.delim(args[2], header = FALSE)
contamination_coords <- data.frame(cont_start = contamination_bed$V2, cont_end = contamination_bed$V3)

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

depcov1 <- ggplot() +
  geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
  labs(title = paste0(id_sample),
       y = "Per base coverage (x)", x = NULL) +
  scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                     expand = expansion(0, 0), limits = c(0, 30000)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  theme_light(base_size = 10) +
  scale_y_log10(
#    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(angle = 90, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(hjust = 1, size = 8)) +
  geom_hline(yintercept = 10, linetype = "dotted", colour = "#5A5A5A")

depcov2 <- ggplot() +
  geom_rect(data = contamination_coords, aes(xmin = cont_start, xmax = cont_end, ymin = 0, ymax = Inf), linewidth = .01, colour = "#5A5A5A", alpha = .1) +
  geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "#000000") +
  labs(title = paste0(id_sample),
       y = "Per base coverage (x)", x = NULL) +
  scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                     expand = expansion(0, 0), limits = c(0, 30000)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  theme_light(base_size = 10) +
  scale_y_log10(
#    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(angle = 90, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(hjust = 1, size = 8)) +
  geom_hline(yintercept = 10, linetype = "dotted", colour = "#5A5A5A")

mapsars1 <- map1 %>% ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
            linewidth = .2, colour = "#000000", alpha = .3) +
  geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
  scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
  theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "UTR" = "#2F67CD",
    "ORFs" = "#FE0B12",
    "Structural proteins" = "#11961B"))

mapsars2 <- map2 %>% ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
            linewidth = .2, colour = "#000000", alpha = .3) +
  geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
  scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
  theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "NSPs" = "#FF8A8E",
    "Accessory factors" = "#D860CF"))

mapsars3 <- map3 %>% ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = class),
            linewidth = .2, colour = "#000000", alpha = .3) +
  geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 3) +
  scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30000)) +
  theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "NSPs" = "#FF8A8E"))

output1 <-  paste0(output, ".coverage.pdf")
plot1 <- depcov1 / mapsars1 / mapsars2 / mapsars3 + plot_layout(nrow = 4, heights = c(3, .4, .3, .3))
save_plot(output1, plot1, base_height = 5, base_width = 16)

output2 <-  paste0(output, ".coverage.contamination.pdf")
plot2 <- depcov2 / mapsars1 / mapsars2 / mapsars3 + plot_layout(nrow = 4, heights = c(3, .4, .3, .3))
save_plot(output2, plot2, base_height = 5, base_width = 16)
