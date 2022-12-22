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

output1 <-  paste0(output, ".coverage.pdf")
save_plot(output1, plot1, base_height = 4, base_width = 16)
