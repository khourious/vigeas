# Program: vigeas [VIral GEnome ASsembly pipelines for WGS]
# Updated: October 02, 2025
# Author: Laise de Moraes <laise.moraes@fiocruz.br>

if(!interactive()) pdf(NULL)

library("ggplot2")
library("cowplot")
library("patchwork")
library("plyr")
library("dplyr")
library("readr")
library("svglite")

args <- commandArgs(TRUE)
id_sample <- args[2]
primer_scheme <- args[3]
output <- args[4]
print(args)

input_depth <- read.delim(args[1], header = FALSE)
depth_coverage <- data.frame(position = input_depth$V2, depth = input_depth$V3)

ref_seq <- list(
  "ARTIC/V1" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ARTIC/V2" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ARTIC/V3" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ARTIC/V4" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ARTIC/V4.1" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ARTIC/V5.3.2" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "FIOCRUZ-IOC/V1" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "FIOCRUZ-IOC/V2" = "Genome reference: MN908947.3 (SARS-CoV-2)",
  "ZikaAsian/V2" = "Genome reference: KJ776791.2 (ZIKV)",
  "DENGUESEQ1/V1" = "Genome reference: NC_001477.1 (DENV-1)",
  "DENGUESEQ2/V1" = "Genome reference: NC_001474.2 (DENV-2)",
  "DENGUESEQ3/V1" = "Genome reference: NC_001475.2 (DENV-3)",
  "DENGUESEQ4/V1" = "Genome reference: NC_002640.1 (DENV-4)",
  "ChikAsianECSA/V1" = "Genome reference: KP164568.1 (CHIKV)",
  "CCEMHTLV1/V1" = "Genome reference: J02029.1 (HTLV-1)",
  "WNV400/V1" = "Genome reference: NC_009942.1 (WNV)",
  "HIV1Sanabani2006/V1" = "Genome reference: K03455.1 (HIV-1)",
  "RSVA/V1" = "Genome reference: MN163126.1 (RSV-A)",
  "RSVA/V2" = "Genome reference: MN163126.1 (RSV-A)",
  "RSVA/V3" = "Genome reference: MN163126.1 (RSV-A)",
  "RSVB/V1" = "Genome reference: MZ515716.1 (RSV-B)",
  "RSVB/V2" = "Genome reference: MZ515716.1 (RSV-B)",
  "RSVB/V3" = "Genome reference: MZ515716.1 (RSV-B)",
  "HTLV1DemincoF/V1" = "Genome reference: J02029.1 (HTLV-1)",
  "DENV1CADDE/V1" = "Genome reference: NC_001477.1 (DENV-1)",
  "DENV2CADDE/V1" = "Genome reference: NC_001474.2 (DENV-2)",
  "DENV2GII2022FN/V1" = "Genome reference: PV789655.1 (DENV-2 GII)",
  "DENV3CADDE/V1" = "Genome reference: NC_001475.2 (DENV-3)",
  "DENV4CADDE/V1" = "Genome reference: NC_002640.1 (DENV-4)",
  "OROVFN400L/V1" = "Genome reference: KP691612.1 (OROV L segment)",
  "OROVFN400M/V1" = "Genome reference: KP691622.1 (OROV M segment)",
  "OROVFN400S/V1" = "Genome reference: KP691623.1 (OROV S segment)")
primer_scheme_2 <- ref_seq[[primer_scheme]]

if (primer_scheme == "ARTIC/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "nCoV-2019_1_LEFT-nCoV-2019_1_RIGHT", 30, 410,
                    "1", "nCoV-2019_3_LEFT-nCoV-2019_3_RIGHT", 642, 1028,
                    "1", "nCoV-2019_5_LEFT-nCoV-2019_5_RIGHT", 1242, 1651,
                    "1", "nCoV-2019_7_LEFT-nCoV-2019_7_RIGHT", 1875, 2269,
                    "1", "nCoV-2019_9_LEFT-nCoV-2019_9_RIGHT", 2505, 2904,
                    "1", "nCoV-2019_11_LEFT-nCoV-2019_11_RIGHT", 3144, 3531,
                    "1", "nCoV-2019_13_LEFT-nCoV-2019_13_RIGHT", 3771, 4164,
                    "1", "nCoV-2019_15_LEFT-nCoV-2019_15_RIGHT", 4294, 4696,
                    "1", "nCoV-2019_17_LEFT-nCoV-2019_17_RIGHT", 4939, 5321,
                    "1", "nCoV-2019_19_LEFT-nCoV-2019_19_RIGHT", 5563, 5957,
                    "1", "nCoV-2019_21_LEFT-nCoV-2019_21_RIGHT", 6167, 6550,
                    "1", "nCoV-2019_23_LEFT-nCoV-2019_23_RIGHT", 6718, 7117,
                    "1", "nCoV-2019_25_LEFT-nCoV-2019_25_RIGHT", 7305, 7694,
                    "1", "nCoV-2019_27_LEFT-nCoV-2019_27_RIGHT", 7943, 8341,
                    "1", "nCoV-2019_29_LEFT-nCoV-2019_29_RIGHT", 8595, 8983,
                    "1", "nCoV-2019_31_LEFT-nCoV-2019_31_RIGHT", 9204, 9585,
                    "1", "nCoV-2019_33_LEFT-nCoV-2019_33_RIGHT", 9784, 10171,
                    "1", "nCoV-2019_35_LEFT-nCoV-2019_35_RIGHT", 10362, 10763,
                    "1", "nCoV-2019_37_LEFT-nCoV-2019_37_RIGHT", 10999, 11394,
                    "1", "nCoV-2019_39_LEFT-nCoV-2019_39_RIGHT", 11555, 11949,
                    "1", "nCoV-2019_41_LEFT-nCoV-2019_41_RIGHT", 12110, 12490,
                    "1", "nCoV-2019_43_LEFT-nCoV-2019_43_RIGHT", 12710, 13096,
                    "1", "nCoV-2019_45_LEFT-nCoV-2019_45_RIGHT", 13319, 13699,
                    "1", "nCoV-2019_47_LEFT-nCoV-2019_47_RIGHT", 13918, 14299,
                    "1", "nCoV-2019_49_LEFT-nCoV-2019_49_RIGHT", 14545, 14926,
                    "1", "nCoV-2019_51_LEFT-nCoV-2019_51_RIGHT", 15171, 15560,
                    "1", "nCoV-2019_53_LEFT-nCoV-2019_53_RIGHT", 15827, 16209,
                    "1", "nCoV-2019_55_LEFT-nCoV-2019_55_RIGHT", 16416, 16833,
                    "1", "nCoV-2019_57_LEFT-nCoV-2019_57_RIGHT", 17065, 17452,
                    "1", "nCoV-2019_59_LEFT-nCoV-2019_59_RIGHT", 17674, 18062,
                    "1", "nCoV-2019_61_LEFT-nCoV-2019_61_RIGHT", 18253, 18672,
                    "1", "nCoV-2019_63_LEFT-nCoV-2019_63_RIGHT", 18896, 19297,
                    "1", "nCoV-2019_65_LEFT-nCoV-2019_65_RIGHT", 19548, 19939,
                    "1", "nCoV-2019_67_LEFT-nCoV-2019_67_RIGHT", 20172, 20572,
                    "1", "nCoV-2019_69_LEFT-nCoV-2019_69_RIGHT", 20786, 21169,
                    "1", "nCoV-2019_71_LEFT-nCoV-2019_71_RIGHT", 21357, 21743,
                    "1", "nCoV-2019_73_LEFT-nCoV-2019_73_RIGHT", 21961, 22346,
                    "1", "nCoV-2019_75_LEFT-nCoV-2019_75_RIGHT", 22516, 22903,
                    "1", "nCoV-2019_77_LEFT-nCoV-2019_77_RIGHT", 23122, 23522,
                    "1", "nCoV-2019_79_LEFT-nCoV-2019_79_RIGHT", 23789, 24169,
                    "1", "nCoV-2019_81_LEFT-nCoV-2019_81_RIGHT", 24391, 24789,
                    "1", "nCoV-2019_83_LEFT-nCoV-2019_83_RIGHT", 24978, 25369,
                    "1", "nCoV-2019_85_LEFT-nCoV-2019_85_RIGHT", 25601, 25994,
                    "1", "nCoV-2019_87_LEFT-nCoV-2019_87_RIGHT", 26197, 26590,
                    "1", "nCoV-2019_89_LEFT-nCoV-2019_89_RIGHT", 26835, 27227,
                    "1", "nCoV-2019_91_LEFT-nCoV-2019_91_RIGHT", 27446, 27854,
                    "1", "nCoV-2019_93_LEFT-nCoV-2019_93_RIGHT", 28081, 28464,
                    "1", "nCoV-2019_95_LEFT-nCoV-2019_95_RIGHT", 28677, 29063,
                    "1", "nCoV-2019_97_LEFT-nCoV-2019_97_RIGHT", 29288, 29693)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "nCoV-2019_2_LEFT-nCoV-2019_2_RIGHT", 320, 726,
                    "2", "nCoV-2019_4_LEFT-nCoV-2019_4_RIGHT", 943, 1337,
                    "2", "nCoV-2019_6_LEFT-nCoV-2019_6_RIGHT", 1573, 1964,
                    "2", "nCoV-2019_8_LEFT-nCoV-2019_8_RIGHT", 2181, 2592,
                    "2", "nCoV-2019_10_LEFT-nCoV-2019_10_RIGHT", 2826, 3210,
                    "2", "nCoV-2019_12_LEFT-nCoV-2019_12_RIGHT", 3460, 3853,
                    "2", "nCoV-2019_14_LEFT-nCoV-2019_14_RIGHT", 4054, 4450,
                    "2", "nCoV-2019_16_LEFT-nCoV-2019_16_RIGHT", 4636, 5017,
                    "2", "nCoV-2019_18_LEFT-nCoV-2019_18_RIGHT", 5230, 5644,
                    "2", "nCoV-2019_20_LEFT-nCoV-2019_20_RIGHT", 5867, 6272,
                    "2", "nCoV-2019_22_LEFT-nCoV-2019_22_RIGHT", 6466, 6873,
                    "2", "nCoV-2019_24_LEFT-nCoV-2019_24_RIGHT", 7035, 7415,
                    "2", "nCoV-2019_26_LEFT-nCoV-2019_26_RIGHT", 7626, 8019,
                    "2", "nCoV-2019_28_LEFT-nCoV-2019_28_RIGHT", 8249, 8661,
                    "2", "nCoV-2019_30_LEFT-nCoV-2019_30_RIGHT", 8888, 9271,
                    "2", "nCoV-2019_32_LEFT-nCoV-2019_32_RIGHT", 9477, 9858,
                    "2", "nCoV-2019_34_LEFT-nCoV-2019_34_RIGHT", 10076, 10459,
                    "2", "nCoV-2019_36_LEFT-nCoV-2019_36_RIGHT", 10666, 11074,
                    "2", "nCoV-2019_38_LEFT-nCoV-2019_38_RIGHT", 11306, 11693,
                    "2", "nCoV-2019_40_LEFT-nCoV-2019_40_RIGHT", 11863, 12256,
                    "2", "nCoV-2019_42_LEFT-nCoV-2019_42_RIGHT", 12417, 12802,
                    "2", "nCoV-2019_44_LEFT-nCoV-2019_44_RIGHT", 13005, 13400,
                    "2", "nCoV-2019_46_LEFT-nCoV-2019_46_RIGHT", 13599, 13984,
                    "2", "nCoV-2019_48_LEFT-nCoV-2019_48_RIGHT", 14207, 14601,
                    "2", "nCoV-2019_50_LEFT-nCoV-2019_50_RIGHT", 14865, 15246,
                    "2", "nCoV-2019_52_LEFT-nCoV-2019_52_RIGHT", 15481, 15886,
                    "2", "nCoV-2019_54_LEFT-nCoV-2019_54_RIGHT", 16118, 16510,
                    "2", "nCoV-2019_56_LEFT-nCoV-2019_56_RIGHT", 16748, 17152,
                    "2", "nCoV-2019_58_LEFT-nCoV-2019_58_RIGHT", 17381, 17761,
                    "2", "nCoV-2019_60_LEFT-nCoV-2019_60_RIGHT", 17966, 18348,
                    "2", "nCoV-2019_62_LEFT-nCoV-2019_62_RIGHT", 18596, 18979,
                    "2", "nCoV-2019_64_LEFT-nCoV-2019_64_RIGHT", 19204, 19616,
                    "2", "nCoV-2019_66_LEFT-nCoV-2019_66_RIGHT", 19844, 20255,
                    "2", "nCoV-2019_68_LEFT-nCoV-2019_68_RIGHT", 20472, 20890,
                    "2", "nCoV-2019_70_LEFT-nCoV-2019_70_RIGHT", 21075, 21455,
                    "2", "nCoV-2019_72_LEFT-nCoV-2019_72_RIGHT", 21658, 22038,
                    "2", "nCoV-2019_74_LEFT-nCoV-2019_74_RIGHT", 22262, 22650,
                    "2", "nCoV-2019_76_LEFT-nCoV-2019_76_RIGHT", 22797, 23214,
                    "2", "nCoV-2019_78_LEFT-nCoV-2019_78_RIGHT", 23443, 23847,
                    "2", "nCoV-2019_80_LEFT-nCoV-2019_80_RIGHT", 24078, 24467,
                    "2", "nCoV-2019_82_LEFT-nCoV-2019_82_RIGHT", 24696, 25076,
                    "2", "nCoV-2019_84_LEFT-nCoV-2019_84_RIGHT", 25279, 25673,
                    "2", "nCoV-2019_86_LEFT-nCoV-2019_86_RIGHT", 25902, 26315,
                    "2", "nCoV-2019_88_LEFT-nCoV-2019_88_RIGHT", 26520, 26913,
                    "2", "nCoV-2019_90_LEFT-nCoV-2019_90_RIGHT", 27141, 27533,
                    "2", "nCoV-2019_92_LEFT-nCoV-2019_92_RIGHT", 27784, 28172,
                    "2", "nCoV-2019_94_LEFT-nCoV-2019_94_RIGHT", 28394, 28779,
                    "2", "nCoV-2019_96_LEFT-nCoV-2019_96_RIGHT", 28985, 29378,
                    "2", "nCoV-2019_98_LEFT-nCoV-2019_98_RIGHT", 29486, 29866)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "FIOCRUZ-IOC/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "hCoV_F1_2kb-hCoV_R1_2kb", 31, 2079,
                    "1", "hCoV_F3_2kb-hCoV_R3_2kb", 3581, 5548,
                    "1", "hCoV_F5_2kb-hCoV_R5_2kb", 7093, 9123,
                    "1", "hCoV_F7_2kb-hCoV_R7_2kb", 10677, 12679,
                    "1", "hCoV_F9_2kb-hCoV_R9_2kb", 17177, 15978,
                    "1", "hCoV_F11_2kb-hCoV_R11_2kb", 17572, 19485,
                    "1", "hCoV_F13_2kb-hCoV_R13_2kb", 21076, 22996,
                    "1", "hCoV_F15_2kb-hCoV_R15_2kb", 24650, 26542,
                    "1", "hCoV_F17_2kb-hCoV_R17_2kb", 27915, 29790)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "hCoV_F2_2kb-hCoV_R2_2kb", 1926, 3737,
                    "2", "hCoV_F4_2kb-hCoV_R4_2kb", 5395, 7255,
                    "2", "hCoV_F6_2kb-hCoV_R6_2kb", 8799, 10837,
                    "2", "hCoV_F8_2kb-hCoV_R8_2kb", 12520, 14328,
                    "2", "hCoV_F10_2kb-hCoV_R10_2kb", 15828, 17754,
                    "2", "hCoV_F12_2kb-hCoV_R12_2kb", 19311, 21241,
                    "2", "hCoV_F14_2kb-hCoV_R14_2kb", 22851, 24812,
                    "2", "hCoV_F16_2kb-hCoV_R16_2kb", 26387, 28351)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ARTIC/V2") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "nCoV-2019_1_LEFT-nCoV-2019_1_RIGHT", 30, 410,
                    "1", "nCoV-2019_3_LEFT-nCoV-2019_3_RIGHT", 642, 1028,
                    "1", "nCoV-2019_5_LEFT-nCoV-2019_5_RIGHT", 1242, 1651,
                    "1", "nCoV-2019_7_LEFT-nCoV-2019_7_RIGHT", 1875, 2269,
                    "1", "nCoV-2019_9_LEFT-nCoV-2019_9_RIGHT", 2505, 2904,
                    "1", "nCoV-2019_11_LEFT-nCoV-2019_11_RIGHT", 3144, 3531,
                    "1", "nCoV-2019_13_LEFT-nCoV-2019_13_RIGHT", 3771, 4164,
                    "1", "nCoV-2019_15_LEFT-nCoV-2019_15_RIGHT", 4294, 4696,
                    "1", "nCoV-2019_17_LEFT-nCoV-2019_17_RIGHT", 4939, 5321,
                    "1", "nCoV-2019_19_LEFT-nCoV-2019_19_RIGHT", 5563, 5957,
                    "1", "nCoV-2019_21_LEFT-nCoV-2019_21_RIGHT", 6167, 6550,
                    "1", "nCoV-2019_23_LEFT-nCoV-2019_23_RIGHT", 6718, 7117,
                    "1", "nCoV-2019_25_LEFT-nCoV-2019_25_RIGHT", 7305, 7694,
                    "1", "nCoV-2019_27_LEFT-nCoV-2019_27_RIGHT", 7943, 8341,
                    "1", "nCoV-2019_29_LEFT-nCoV-2019_29_RIGHT", 8595, 8983,
                    "1", "nCoV-2019_31_LEFT-nCoV-2019_31_RIGHT", 9204, 9585,
                    "1", "nCoV-2019_33_LEFT-nCoV-2019_33_RIGHT", 9784, 10171,
                    "1", "nCoV-2019_35_LEFT-nCoV-2019_35_RIGHT", 10362, 10763,
                    "1", "nCoV-2019_37_LEFT-nCoV-2019_37_RIGHT", 10999, 11394,
                    "1", "nCoV-2019_39_LEFT-nCoV-2019_39_RIGHT", 11555, 11949,
                    "1", "nCoV-2019_41_LEFT-nCoV-2019_41_RIGHT", 12110, 12490,
                    "1", "nCoV-2019_43_LEFT-nCoV-2019_43_RIGHT", 12710, 13096,
                    "1", "nCoV-2019_45_LEFT-nCoV-2019_45_RIGHT", 13319, 13699,
                    "1", "nCoV-2019_47_LEFT-nCoV-2019_47_RIGHT", 13918, 14299,
                    "1", "nCoV-2019_49_LEFT-nCoV-2019_49_RIGHT", 14545, 14926,
                    "1", "nCoV-2019_51_LEFT-nCoV-2019_51_RIGHT", 15171, 15560,
                    "1", "nCoV-2019_53_LEFT-nCoV-2019_53_RIGHT", 15827, 16209,
                    "1", "nCoV-2019_55_LEFT-nCoV-2019_55_RIGHT", 16416, 16833,
                    "1", "nCoV-2019_57_LEFT-nCoV-2019_57_RIGHT", 17065, 17452,
                    "1", "nCoV-2019_59_LEFT-nCoV-2019_59_RIGHT", 17674, 18062,
                    "1", "nCoV-2019_61_LEFT-nCoV-2019_61_RIGHT", 18253, 18672,
                    "1", "nCoV-2019_63_LEFT-nCoV-2019_63_RIGHT", 18896, 19297,
                    "1", "nCoV-2019_65_LEFT-nCoV-2019_65_RIGHT", 19548, 19939,
                    "1", "nCoV-2019_67_LEFT-nCoV-2019_67_RIGHT", 20172, 20572,
                    "1", "nCoV-2019_69_LEFT-nCoV-2019_69_RIGHT", 20786, 21169,
                    "1", "nCoV-2019_71_LEFT-nCoV-2019_71_RIGHT", 21357, 21743,
                    "1", "nCoV-2019_73_LEFT-nCoV-2019_73_RIGHT", 21961, 22346,
                    "1", "nCoV-2019_75_LEFT-nCoV-2019_75_RIGHT", 22516, 22903,
                    "1", "nCoV-2019_77_LEFT-nCoV-2019_77_RIGHT", 23122, 23522,
                    "1", "nCoV-2019_79_LEFT-nCoV-2019_79_RIGHT", 23789, 24169,
                    "1", "nCoV-2019_81_LEFT-nCoV-2019_81_RIGHT", 24391, 24789,
                    "1", "nCoV-2019_83_LEFT-nCoV-2019_83_RIGHT", 24978, 25369,
                    "1", "nCoV-2019_85_LEFT-nCoV-2019_85_RIGHT", 25601, 25994,
                    "1", "nCoV-2019_87_LEFT-nCoV-2019_87_RIGHT", 26197, 26590,
                    "1", "nCoV-2019_89_LEFT-nCoV-2019_89_RIGHT", 26835, 27227,
                    "1", "nCoV-2019_91_LEFT-nCoV-2019_91_RIGHT", 27446, 27854,
                    "1", "nCoV-2019_93_LEFT-nCoV-2019_93_RIGHT", 28081, 28464,
                    "1", "nCoV-2019_95_LEFT-nCoV-2019_95_RIGHT", 28677, 29063,
                    "1", "nCoV-2019_97_LEFT-nCoV-2019_97_RIGHT", 29288, 29693)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "nCoV-2019_2_LEFT-nCoV-2019_2_RIGHT", 320, 726,
                    "2", "nCoV-2019_4_LEFT-nCoV-2019_4_RIGHT", 943, 1337,
                    "2", "nCoV-2019_6_LEFT-nCoV-2019_6_RIGHT", 1573, 1964,
                    "2", "nCoV-2019_8_LEFT-nCoV-2019_8_RIGHT", 2181, 2592,
                    "2", "nCoV-2019_10_LEFT-nCoV-2019_10_RIGHT", 2826, 3210,
                    "2", "nCoV-2019_12_LEFT-nCoV-2019_12_RIGHT", 3460, 3853,
                    "2", "nCoV-2019_14_LEFT-nCoV-2019_14_RIGHT", 4054, 4450,
                    "2", "nCoV-2019_16_LEFT-nCoV-2019_16_RIGHT", 4636, 5017,
                    "2", "nCoV-2019_18_LEFT_alt2-nCoV-2019_18_RIGHT", 5257, 5644,
                    "2", "nCoV-2019_20_LEFT-nCoV-2019_20_RIGHT", 5867, 6272,
                    "2", "nCoV-2019_22_LEFT-nCoV-2019_22_RIGHT", 6466, 6873,
                    "2", "nCoV-2019_24_LEFT-nCoV-2019_24_RIGHT", 7035, 7415,
                    "2", "nCoV-2019_26_LEFT-nCoV-2019_26_RIGHT", 7626, 8019,
                    "2", "nCoV-2019_28_LEFT-nCoV-2019_28_RIGHT", 8249, 8661,
                    "2", "nCoV-2019_30_LEFT-nCoV-2019_30_RIGHT", 8888, 9271,
                    "2", "nCoV-2019_32_LEFT-nCoV-2019_32_RIGHT", 9477, 9858,
                    "2", "nCoV-2019_34_LEFT-nCoV-2019_34_RIGHT", 10076, 10459,
                    "2", "nCoV-2019_36_LEFT-nCoV-2019_36_RIGHT", 10666, 11074,
                    "2", "nCoV-2019_38_LEFT-nCoV-2019_38_RIGHT", 11306, 11693,
                    "2", "nCoV-2019_40_LEFT-nCoV-2019_40_RIGHT", 11863, 12256,
                    "2", "nCoV-2019_42_LEFT-nCoV-2019_42_RIGHT", 12417, 12802,
                    "2", "nCoV-2019_44_LEFT-nCoV-2019_44_RIGHT", 13005, 13400,
                    "2", "nCoV-2019_46_LEFT-nCoV-2019_46_RIGHT", 13599, 13984,
                    "2", "nCoV-2019_48_LEFT-nCoV-2019_48_RIGHT", 14207, 14601,
                    "2", "nCoV-2019_50_LEFT-nCoV-2019_50_RIGHT", 14865, 15246,
                    "2", "nCoV-2019_52_LEFT-nCoV-2019_52_RIGHT", 15481, 15886,
                    "2", "nCoV-2019_54_LEFT-nCoV-2019_54_RIGHT", 16118, 16510,
                    "2", "nCoV-2019_56_LEFT-nCoV-2019_56_RIGHT", 16748, 17152,
                    "2", "nCoV-2019_58_LEFT-nCoV-2019_58_RIGHT", 17381, 17761,
                    "2", "nCoV-2019_60_LEFT-nCoV-2019_60_RIGHT", 17966, 18348,
                    "2", "nCoV-2019_62_LEFT-nCoV-2019_62_RIGHT", 18596, 18979,
                    "2", "nCoV-2019_64_LEFT-nCoV-2019_64_RIGHT", 19204, 19616,
                    "2", "nCoV-2019_66_LEFT-nCoV-2019_66_RIGHT", 19844, 20255,
                    "2", "nCoV-2019_68_LEFT-nCoV-2019_68_RIGHT", 20472, 20890,
                    "2", "nCoV-2019_70_LEFT-nCoV-2019_70_RIGHT", 21075, 21455,
                    "2", "nCoV-2019_72_LEFT-nCoV-2019_72_RIGHT", 21658, 22038,
                    "2", "nCoV-2019_74_LEFT-nCoV-2019_74_RIGHT", 22262, 22650,
                    "2", "nCoV-2019_76_LEFT-nCoV-2019_76_RIGHT", 22797, 23214,
                    "2", "nCoV-2019_78_LEFT-nCoV-2019_78_RIGHT", 23443, 23847,
                    "2", "nCoV-2019_80_LEFT-nCoV-2019_80_RIGHT", 24078, 24467,
                    "2", "nCoV-2019_82_LEFT-nCoV-2019_82_RIGHT", 24696, 25076,
                    "2", "nCoV-2019_84_LEFT-nCoV-2019_84_RIGHT", 25279, 25673,
                    "2", "nCoV-2019_86_LEFT-nCoV-2019_86_RIGHT", 25902, 26315,
                    "2", "nCoV-2019_88_LEFT-nCoV-2019_88_RIGHT", 26520, 26913,
                    "2", "nCoV-2019_90_LEFT-nCoV-2019_90_RIGHT", 27141, 27533,
                    "2", "nCoV-2019_92_LEFT-nCoV-2019_92_RIGHT", 27784, 28172,
                    "2", "nCoV-2019_94_LEFT-nCoV-2019_94_RIGHT", 28394, 28779,
                    "2", "nCoV-2019_96_LEFT-nCoV-2019_96_RIGHT", 28985, 29378,
                    "2", "nCoV-2019_98_LEFT-nCoV-2019_98_RIGHT", 29486, 29866)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "FIOCRUZ-IOC/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "hCoV_F1_2kb-hCoV_R1_2kb", 31, 2079,
                    "1", "hCoV_F3_2kb-hCoV_R3_2kb", 3581, 5548,
                    "1", "hCoV_F5_2kb-hCoV_R5_2kb", 7093, 9123,
                    "1", "hCoV_F7_2kb-hCoV_R7_2kb", 10677, 12679,
                    "1", "hCoV_F9_2kb-hCoV_R9_2kb", 17177, 15978,
                    "1", "hCoV_F11_2kb-hCoV_R11_2kb", 17572, 19485,
                    "1", "hCoV_F13_2kb-hCoV_R13_2kb", 21076, 22996,
                    "1", "hCoV_F15_2kb-hCoV_R15_2kb", 24650, 26542,
                    "1", "hCoV_F17_2kb-hCoV_R17_2kb", 27915, 29790)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "hCoV_F2_2kb-hCoV_R2_2kb", 1926, 3737,
                    "2", "hCoV_F4_2kb-hCoV_R4_2kb", 5395, 7255,
                    "2", "hCoV_F6_2kb-hCoV_R6_2kb", 8799, 10837,
                    "2", "hCoV_F8_2kb-hCoV_R8_2kb", 12520, 14328,
                    "2", "hCoV_F10_2kb-hCoV_R10_2kb", 15828, 17754,
                    "2", "hCoV_F12_2kb-hCoV_R12_2kb", 19311, 21241,
                    "2", "hCoV_F14_2kb-hCoV_R14_2kb", 22851, 24812,
                    "2", "hCoV_F16_2kb-hCoV_R16_2kb", 26387, 28351)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                         # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ARTIC/V3") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "nCoV-2019_1_LEFT-nCoV-2019_1_RIGHT", 30, 410,
                    "1", "nCoV-2019_3_LEFT-nCoV-2019_3_RIGHT", 642, 1028,
                    "1", "nCoV-2019_5_LEFT-nCoV-2019_5_RIGHT", 1242, 1651,
                    "1", "nCoV-2019_7_LEFT-nCoV-2019_7_RIGHT", 1875, 2269,
                    "1", "nCoV-2019_9_LEFT-nCoV-2019_9_RIGHT", 2505, 2904,
                    "1", "nCoV-2019_11_LEFT-nCoV-2019_11_RIGHT", 3144, 3531,
                    "1", "nCoV-2019_13_LEFT-nCoV-2019_13_RIGHT", 3771, 4164,
                    "1", "nCoV-2019_15_LEFT-nCoV-2019_15_RIGHT", 4294, 4696,
                    "1", "nCoV-2019_17_LEFT-nCoV-2019_17_RIGHT", 4939, 5321,
                    "1", "nCoV-2019_19_LEFT-nCoV-2019_19_RIGHT", 5563, 5957,
                    "1", "nCoV-2019_21_LEFT-nCoV-2019_21_RIGHT", 6167, 6550,
                    "1", "nCoV-2019_23_LEFT-nCoV-2019_23_RIGHT", 6718, 7117,
                    "1", "nCoV-2019_25_LEFT-nCoV-2019_25_RIGHT", 7305, 7694,
                    "1", "nCoV-2019_27_LEFT-nCoV-2019_27_RIGHT", 7943, 8341,
                    "1", "nCoV-2019_29_LEFT-nCoV-2019_29_RIGHT", 8595, 8983,
                    "1", "nCoV-2019_31_LEFT-nCoV-2019_31_RIGHT", 9204, 9585,
                    "1", "nCoV-2019_33_LEFT-nCoV-2019_33_RIGHT", 9784, 10171,
                    "1", "nCoV-2019_35_LEFT-nCoV-2019_35_RIGHT", 10362, 10763,
                    "1", "nCoV-2019_37_LEFT-nCoV-2019_37_RIGHT", 10999, 11394,
                    "1", "nCoV-2019_39_LEFT-nCoV-2019_39_RIGHT", 11555, 11949,
                    "1", "nCoV-2019_41_LEFT-nCoV-2019_41_RIGHT", 12110, 12490,
                    "1", "nCoV-2019_43_LEFT-nCoV-2019_43_RIGHT", 12710, 13096,
                    "1", "nCoV-2019_45_LEFT-nCoV-2019_45_RIGHT", 13319, 13699,
                    "1", "nCoV-2019_47_LEFT-nCoV-2019_47_RIGHT", 13918, 14299,
                    "1", "nCoV-2019_49_LEFT-nCoV-2019_49_RIGHT", 14545, 14926,
                    "1", "nCoV-2019_51_LEFT-nCoV-2019_51_RIGHT", 15171, 15560,
                    "1", "nCoV-2019_53_LEFT-nCoV-2019_53_RIGHT", 15827, 16209,
                    "1", "nCoV-2019_55_LEFT-nCoV-2019_55_RIGHT", 16416, 16833,
                    "1", "nCoV-2019_57_LEFT-nCoV-2019_57_RIGHT", 17065, 17452,
                    "1", "nCoV-2019_59_LEFT-nCoV-2019_59_RIGHT", 17674, 18062,
                    "1", "nCoV-2019_61_LEFT-nCoV-2019_61_RIGHT", 18253, 18672,
                    "1", "nCoV-2019_63_LEFT-nCoV-2019_63_RIGHT", 18896, 19297,
                    "1", "nCoV-2019_65_LEFT-nCoV-2019_65_RIGHT", 19548, 19939,
                    "1", "nCoV-2019_67_LEFT-nCoV-2019_67_RIGHT", 20172, 20572,
                    "1", "nCoV-2019_69_LEFT-nCoV-2019_69_RIGHT", 20786, 21169,
                    "1", "nCoV-2019_71_LEFT-nCoV-2019_71_RIGHT", 21357, 21743,
                    "1", "nCoV-2019_73_LEFT-nCoV-2019_73_RIGHT", 21961, 22346,
                    "1", "nCoV-2019_75_LEFT-nCoV-2019_75_RIGHT", 22516, 22903,
                    "1", "nCoV-2019_77_LEFT-nCoV-2019_77_RIGHT", 23122, 23522,
                    "1", "nCoV-2019_79_LEFT-nCoV-2019_79_RIGHT", 23789, 24169,
                    "1", "nCoV-2019_81_LEFT-nCoV-2019_81_RIGHT", 24391, 24789,
                    "1", "nCoV-2019_83_LEFT-nCoV-2019_83_RIGHT", 24978, 25369,
                    "1", "nCoV-2019_85_LEFT-nCoV-2019_85_RIGHT", 25601, 25994,
                    "1", "nCoV-2019_87_LEFT-nCoV-2019_87_RIGHT", 26197, 26590,
                    "1", "nCoV-2019_89_LEFT-nCoV-2019_89_RIGHT", 26835, 27227,
                    "1", "nCoV-2019_91_LEFT-nCoV-2019_91_RIGHT", 27446, 27854,
                    "1", "nCoV-2019_93_LEFT-nCoV-2019_93_RIGHT", 28081, 28464,
                    "1", "nCoV-2019_95_LEFT-nCoV-2019_95_RIGHT", 28677, 29063,
                    "1", "nCoV-2019_97_LEFT-nCoV-2019_97_RIGHT", 29288, 29693)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "nCoV-2019_7_LEFT_alt0-nCoV-2019_7_RIGHT_alt5", 1868, 2264,
                    "1", "nCoV-2019_9_LEFT_alt4-nCoV-2019_9_RIGHT_alt2", 2504, 2902,
                    "1", "nCoV-2019_15_LEFT_alt1-nCoV-2019_15_RIGHT_alt3", 4296, 4689,
                    "1", "nCoV-2019_21_LEFT_alt2-nCoV-2019_21_RIGHT_alt0", 6168, 6548,
                    "1", "nCoV-2019_45_LEFT_alt2-nCoV-2019_45_RIGHT_alt7", 13307, 13689,
                    "1", "nCoV-2019_89_LEFT_alt2-nCoV-2019_89_RIGHT_alt4", 26838, 27215)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p3 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "nCoV-2019_2_LEFT-nCoV-2019_2_RIGHT", 320, 726,
                    "2", "nCoV-2019_4_LEFT-nCoV-2019_4_RIGHT", 943, 1337,
                    "2", "nCoV-2019_6_LEFT-nCoV-2019_6_RIGHT", 1573, 1964,
                    "2", "nCoV-2019_8_LEFT-nCoV-2019_8_RIGHT", 2181, 2592,
                    "2", "nCoV-2019_10_LEFT-nCoV-2019_10_RIGHT", 2826, 3210,
                    "2", "nCoV-2019_12_LEFT-nCoV-2019_12_RIGHT", 3460, 3853,
                    "2", "nCoV-2019_14_LEFT-nCoV-2019_14_RIGHT", 4054, 4450,
                    "2", "nCoV-2019_16_LEFT-nCoV-2019_16_RIGHT", 4636, 5017,
                    "2", "nCoV-2019_18_LEFT-nCoV-2019_18_RIGHT", 5230, 5644,
                    "2", "nCoV-2019_20_LEFT-nCoV-2019_20_RIGHT", 5867, 6272,
                    "2", "nCoV-2019_22_LEFT-nCoV-2019_22_RIGHT", 6466, 6873,
                    "2", "nCoV-2019_24_LEFT-nCoV-2019_24_RIGHT", 7035, 7415,
                    "2", "nCoV-2019_26_LEFT-nCoV-2019_26_RIGHT", 7626, 8019,
                    "2", "nCoV-2019_28_LEFT-nCoV-2019_28_RIGHT", 8249, 8661,
                    "2", "nCoV-2019_30_LEFT-nCoV-2019_30_RIGHT", 8888, 9271,
                    "2", "nCoV-2019_32_LEFT-nCoV-2019_32_RIGHT", 9477, 9858,
                    "2", "nCoV-2019_34_LEFT-nCoV-2019_34_RIGHT", 10076, 10459,
                    "2", "nCoV-2019_36_LEFT-nCoV-2019_36_RIGHT", 10666, 11074,
                    "2", "nCoV-2019_38_LEFT-nCoV-2019_38_RIGHT", 11306, 11693,
                    "2", "nCoV-2019_40_LEFT-nCoV-2019_40_RIGHT", 11863, 12256,
                    "2", "nCoV-2019_42_LEFT-nCoV-2019_42_RIGHT", 12417, 12802,
                    "2", "nCoV-2019_44_LEFT-nCoV-2019_44_RIGHT", 13005, 13400,
                    "2", "nCoV-2019_46_LEFT-nCoV-2019_46_RIGHT", 13599, 13984,
                    "2", "nCoV-2019_48_LEFT-nCoV-2019_48_RIGHT", 14207, 14601,
                    "2", "nCoV-2019_50_LEFT-nCoV-2019_50_RIGHT", 14865, 15246,
                    "2", "nCoV-2019_52_LEFT-nCoV-2019_52_RIGHT", 15481, 15886,
                    "2", "nCoV-2019_54_LEFT-nCoV-2019_54_RIGHT", 16118, 16510,
                    "2", "nCoV-2019_56_LEFT-nCoV-2019_56_RIGHT", 16748, 17152,
                    "2", "nCoV-2019_58_LEFT-nCoV-2019_58_RIGHT", 17381, 17761,
                    "2", "nCoV-2019_60_LEFT-nCoV-2019_60_RIGHT", 17966, 18348,
                    "2", "nCoV-2019_62_LEFT-nCoV-2019_62_RIGHT", 18596, 18979,
                    "2", "nCoV-2019_64_LEFT-nCoV-2019_64_RIGHT", 19204, 19616,
                    "2", "nCoV-2019_66_LEFT-nCoV-2019_66_RIGHT", 19844, 20255,
                    "2", "nCoV-2019_68_LEFT-nCoV-2019_68_RIGHT", 20472, 20890,
                    "2", "nCoV-2019_70_LEFT-nCoV-2019_70_RIGHT", 21075, 21455,
                    "2", "nCoV-2019_72_LEFT-nCoV-2019_72_RIGHT", 21658, 22038,
                    "2", "nCoV-2019_74_LEFT-nCoV-2019_74_RIGHT", 22262, 22650,
                    "2", "nCoV-2019_76_LEFT-nCoV-2019_76_RIGHT", 22797, 23214,
                    "2", "nCoV-2019_78_LEFT-nCoV-2019_78_RIGHT", 23443, 23847,
                    "2", "nCoV-2019_80_LEFT-nCoV-2019_80_RIGHT", 24078, 24467,
                    "2", "nCoV-2019_82_LEFT-nCoV-2019_82_RIGHT", 24696, 25076,
                    "2", "nCoV-2019_84_LEFT-nCoV-2019_84_RIGHT", 25279, 25673,
                    "2", "nCoV-2019_86_LEFT-nCoV-2019_86_RIGHT", 25902, 26315,
                    "2", "nCoV-2019_88_LEFT-nCoV-2019_88_RIGHT", 26520, 26913,
                    "2", "nCoV-2019_90_LEFT-nCoV-2019_90_RIGHT", 27141, 27533,
                    "2", "nCoV-2019_92_LEFT-nCoV-2019_92_RIGHT", 27784, 28172,
                    "2", "nCoV-2019_94_LEFT-nCoV-2019_94_RIGHT", 28394, 28779,
                    "2", "nCoV-2019_96_LEFT-nCoV-2019_96_RIGHT", 28985, 29378,
                    "2", "nCoV-2019_98_LEFT-nCoV-2019_98_RIGHT", 29486, 29866)
  map1plot3 <- map1p3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map1p4 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "nCoV-2019_14_LEFT_alt4-nCoV-2019_14_RIGHT_alt2", 4044, 4424,
                    "2", "nCoV-2019_18_LEFT_alt2-nCoV-2019_18_RIGHT_alt1", 5257, 5643,
                    "2", "nCoV-2019_44_LEFT_alt3-nCoV-2019_44_RIGHT_alt0", 13007, 13385,
                    "2", "nCoV-2019_46_LEFT_alt1-nCoV-2019_46_RIGHT_alt2", 13602, 13984,
                    "2", "nCoV-2019_76_LEFT_alt3-nCoV-2019_76_RIGHT_alt0", 22798, 23212)
  map1plot4 <- map1p4 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / map1plot3 / map1plot4 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 9, heights = c(3, .1, .1, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ARTIC/V4") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "SARS-CoV-2_1_LEFT-SARS-CoV-2_1_RIGHT", 25, 431,
                    "1", "SARS-CoV-2_3_LEFT-SARS-CoV-2_3_RIGHT", 644, 1044,
                    "1", "SARS-CoV-2_5_LEFT-SARS-CoV-2_5_RIGHT", 1245, 1650,
                    "1", "SARS-CoV-2_7_LEFT-SARS-CoV-2_7_RIGHT", 1851, 2250,
                    "1", "SARS-CoV-2_9_LEFT-SARS-CoV-2_9_RIGHT", 2483, 2885,
                    "1", "SARS-CoV-2_11_LEFT-SARS-CoV-2_11_RIGHT", 3078, 3492,
                    "1", "SARS-CoV-2_13_LEFT-SARS-CoV-2_13_RIGHT", 3683, 4093,
                    "1", "SARS-CoV-2_15_LEFT-SARS-CoV-2_15_RIGHT", 4312, 4710,
                    "1", "SARS-CoV-2_17_LEFT-SARS-CoV-2_17_RIGHT", 4923, 5331,
                    "1", "SARS-CoV-2_19_LEFT-SARS-CoV-2_19_RIGHT", 5561, 5957,
                    "1", "SARS-CoV-2_21_LEFT-SARS-CoV-2_21_RIGHT", 6184, 6582,
                    "1", "SARS-CoV-2_23_LEFT-SARS-CoV-2_23_RIGHT", 6747, 7148,
                    "1", "SARS-CoV-2_25_LEFT-SARS-CoV-2_25_RIGHT", 7381, 7770,
                    "1", "SARS-CoV-2_27_LEFT-SARS-CoV-2_27_RIGHT", 7997, 8395,
                    "1", "SARS-CoV-2_29_LEFT-SARS-CoV-2_29_RIGHT", 8596, 9013,
                    "1", "SARS-CoV-2_31_LEFT-SARS-CoV-2_31_RIGHT", 9168, 9564,
                    "1", "SARS-CoV-2_33_LEFT-SARS-CoV-2_33_RIGHT", 9782, 10176,
                    "1", "SARS-CoV-2_35_LEFT-SARS-CoV-2_35_RIGHT", 10393, 10810,
                    "1", "SARS-CoV-2_37_LEFT-SARS-CoV-2_37_RIGHT", 11000, 11414,
                    "1", "SARS-CoV-2_39_LEFT-SARS-CoV-2_39_RIGHT", 11624, 12033,
                    "1", "SARS-CoV-2_41_LEFT-SARS-CoV-2_41_RIGHT", 12234, 12643,
                    "1", "SARS-CoV-2_43_LEFT-SARS-CoV-2_43_RIGHT", 12831, 13240,
                    "1", "SARS-CoV-2_45_LEFT-SARS-CoV-2_45_RIGHT", 13463, 13859,
                    "1", "SARS-CoV-2_47_LEFT-SARS-CoV-2_47_RIGHT", 14045, 14457,
                    "1", "SARS-CoV-2_49_LEFT-SARS-CoV-2_49_RIGHT", 14647, 15050,
                    "1", "SARS-CoV-2_51_LEFT-SARS-CoV-2_51_RIGHT", 15214, 15619,
                    "1", "SARS-CoV-2_53_LEFT-SARS-CoV-2_53_RIGHT", 15855, 16260,
                    "1", "SARS-CoV-2_55_LEFT-SARS-CoV-2_55_RIGHT", 16386, 16796,
                    "1", "SARS-CoV-2_57_LEFT-SARS-CoV-2_57_RIGHT", 16986, 17405,
                    "1", "SARS-CoV-2_59_LEFT-SARS-CoV-2_59_RIGHT", 17615, 18022,
                    "1", "SARS-CoV-2_61_LEFT-SARS-CoV-2_61_RIGHT", 18244, 18652,
                    "1", "SARS-CoV-2_63_LEFT-SARS-CoV-2_63_RIGHT", 18869, 19277,
                    "1", "SARS-CoV-2_65_LEFT-SARS-CoV-2_65_RIGHT", 19485, 19901,
                    "1", "SARS-CoV-2_67_LEFT-SARS-CoV-2_67_RIGHT", 20090, 20497,
                    "1", "SARS-CoV-2_69_LEFT-SARS-CoV-2_69_RIGHT", 20677, 21080,
                    "1", "SARS-CoV-2_71_LEFT-SARS-CoV-2_71_RIGHT", 21294, 21700,
                    "1", "SARS-CoV-2_73_LEFT-SARS-CoV-2_73_RIGHT", 21865, 22274,
                    "1", "SARS-CoV-2_75_LEFT-SARS-CoV-2_75_RIGHT", 22402, 22805,
                    "1", "SARS-CoV-2_77_LEFT-SARS-CoV-2_77_RIGHT", 22944, 23351,
                    "1", "SARS-CoV-2_79_LEFT-SARS-CoV-2_79_RIGHT", 23553, 23955,
                    "1", "SARS-CoV-2_81_LEFT-SARS-CoV-2_81_RIGHT", 24171, 24567,
                    "1", "SARS-CoV-2_83_LEFT-SARS-CoV-2_83_RIGHT", 24750, 25150,
                    "1", "SARS-CoV-2_85_LEFT-SARS-CoV-2_85_RIGHT", 25331, 25740,
                    "1", "SARS-CoV-2_87_LEFT-SARS-CoV-2_87_RIGHT", 25951, 26360,
                    "1", "SARS-CoV-2_89_LEFT-SARS-CoV-2_89_RIGHT", 26564, 26979,
                    "1", "SARS-CoV-2_91_LEFT-SARS-CoV-2_91_RIGHT", 27152, 27560,
                    "1", "SARS-CoV-2_93_LEFT-SARS-CoV-2_93_RIGHT", 27700, 28104,
                    "1", "SARS-CoV-2_95_LEFT-SARS-CoV-2_95_RIGHT", 28190, 28598,
                    "1", "SARS-CoV-2_97_LEFT-SARS-CoV-2_97_RIGHT", 28827, 29227,
                    "1", "SARS-CoV-2_99_LEFT-SARS-CoV-2_99_RIGHT", 29452, 29854)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "SARS-CoV-2_2_LEFT-SARS-CoV-2_2_RIGHT", 324, 727,
                    "2", "SARS-CoV-2_4_LEFT-SARS-CoV-2_4_RIGHT", 944, 1362,
                    "2", "SARS-CoV-2_6_LEFT-SARS-CoV-2_6_RIGHT", 1540, 1948,
                    "2", "SARS-CoV-2_8_LEFT-SARS-CoV-2_8_RIGHT", 2154, 2571,
                    "2", "SARS-CoV-2_10_LEFT-SARS-CoV-2_10_RIGHT", 2826, 3210,
                    "2", "SARS-CoV-2_12_LEFT-SARS-CoV-2_12_RIGHT", 3390, 3794,
                    "2", "SARS-CoV-2_14_LEFT-SARS-CoV-2_14_RIGHT", 3992, 4409,
                    "2", "SARS-CoV-2_16_LEFT-SARS-CoV-2_16_RIGHT", 4620, 5017,
                    "2", "SARS-CoV-2_18_LEFT-SARS-CoV-2_18_RIGHT", 5230, 5643,
                    "2", "SARS-CoV-2_20_LEFT-SARS-CoV-2_20_RIGHT", 5867, 6272,
                    "2", "SARS-CoV-2_22_LEFT-SARS-CoV-2_22_RIGHT", 6478, 6885,
                    "2", "SARS-CoV-2_24_LEFT-SARS-CoV-2_24_RIGHT", 7057, 7467,
                    "2", "SARS-CoV-2_26_LEFT-SARS-CoV-2_26_RIGHT", 7672, 8092,
                    "2", "SARS-CoV-2_28_LEFT-SARS-CoV-2_28_RIGHT", 8304, 8714,
                    "2", "SARS-CoV-2_30_LEFT-SARS-CoV-2_30_RIGHT", 8919, 9329,
                    "2", "SARS-CoV-2_32_LEFT-SARS-CoV-2_32_RIGHT", 9470, 9866,
                    "2", "SARS-CoV-2_34_LEFT-SARS-CoV-2_34_RIGHT", 10076, 10491,
                    "2", "SARS-CoV-2_36_LEFT-SARS-CoV-2_36_RIGHT", 10713, 11116,
                    "2", "SARS-CoV-2_38_LEFT-SARS-CoV-2_38_RIGHT", 11305, 11720,
                    "2", "SARS-CoV-2_40_LEFT-SARS-CoV-2_40_RIGHT", 11937, 12339,
                    "2", "SARS-CoV-2_42_LEFT-SARS-CoV-2_42_RIGHT", 12519, 12920,
                    "2", "SARS-CoV-2_44_LEFT-SARS-CoV-2_44_RIGHT", 13124, 13528,
                    "2", "SARS-CoV-2_46_LEFT-SARS-CoV-2_46_RIGHT", 13752, 14144,
                    "2", "SARS-CoV-2_48_LEFT-SARS-CoV-2_48_RIGHT", 14338, 14743,
                    "2", "SARS-CoV-2_50_LEFT-SARS-CoV-2_50_RIGHT", 14953, 15358,
                    "2", "SARS-CoV-2_52_LEFT-SARS-CoV-2_52_RIGHT", 15535, 15941,
                    "2", "SARS-CoV-2_54_LEFT-SARS-CoV-2_54_RIGHT", 16112, 16508,
                    "2", "SARS-CoV-2_56_LEFT-SARS-CoV-2_56_RIGHT", 16692, 17105,
                    "2", "SARS-CoV-2_58_LEFT-SARS-CoV-2_58_RIGHT", 17323, 17711,
                    "2", "SARS-CoV-2_60_LEFT-SARS-CoV-2_60_RIGHT", 17911, 18328,
                    "2", "SARS-CoV-2_62_LEFT-SARS-CoV-2_62_RIGHT", 18550, 18961,
                    "2", "SARS-CoV-2_64_LEFT-SARS-CoV-2_64_RIGHT", 19183, 19586,
                    "2", "SARS-CoV-2_66_LEFT-SARS-CoV-2_66_RIGHT", 19810, 20216,
                    "2", "SARS-CoV-2_68_LEFT-SARS-CoV-2_68_RIGHT", 20377, 20792,
                    "2", "SARS-CoV-2_70_LEFT-SARS-CoV-2_70_RIGHT", 20988, 21387,
                    "2", "SARS-CoV-2_72_LEFT-SARS-CoV-2_72_RIGHT", 21532, 21933,
                    "2", "SARS-CoV-2_74_LEFT-SARS-CoV-2_74_RIGHT", 22091, 22503,
                    "2", "SARS-CoV-2_76_LEFT-SARS-CoV-2_76_RIGHT", 22648, 23057,
                    "2", "SARS-CoV-2_78_LEFT-SARS-CoV-2_78_RIGHT", 23219, 23635,
                    "2", "SARS-CoV-2_80_LEFT-SARS-CoV-2_80_RIGHT", 23853, 24258,
                    "2", "SARS-CoV-2_82_LEFT-SARS-CoV-2_82_RIGHT", 24426, 24836,
                    "2", "SARS-CoV-2_84_LEFT-SARS-CoV-2_84_RIGHT", 25051, 25461,
                    "2", "SARS-CoV-2_86_LEFT-SARS-CoV-2_86_RIGHT", 25645, 26050,
                    "2", "SARS-CoV-2_88_LEFT-SARS-CoV-2_88_RIGHT", 26255, 26661,
                    "2", "SARS-CoV-2_90_LEFT-SARS-CoV-2_90_RIGHT", 26873, 27283,
                    "2", "SARS-CoV-2_92_LEFT-SARS-CoV-2_92_RIGHT", 27447, 27855,
                    "2", "SARS-CoV-2_94_LEFT-SARS-CoV-2_94_RIGHT", 27996, 28416,
                    "2", "SARS-CoV-2_96_LEFT-SARS-CoV-2_96_RIGHT", 28512, 28914,
                    "2", "SARS-CoV-2_98_LEFT-SARS-CoV-2_98_RIGHT", 29136, 29534)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ARTIC/V4.1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "SARS-CoV-2_1_LEFT-SARS-CoV-2_1_RIGHT", 25, 431,
                    "1", "SARS-CoV-2_3_LEFT-SARS-CoV-2_3_RIGHT", 644, 1044,
                    "1", "SARS-CoV-2_5_LEFT-SARS-CoV-2_5_RIGHT", 1245, 1650,
                    "1", "SARS-CoV-2_7_LEFT-SARS-CoV-2_7_RIGHT", 1851, 2250,
                    "1", "SARS-CoV-2_9_LEFT-SARS-CoV-2_9_RIGHT", 2483, 2885,
                    "1", "SARS-CoV-2_11_LEFT-SARS-CoV-2_11_RIGHT", 3078, 3492,
                    "1", "SARS-CoV-2_13_LEFT-SARS-CoV-2_13_RIGHT", 3683, 4093,
                    "1", "SARS-CoV-2_15_LEFT-SARS-CoV-2_15_RIGHT", 4312, 4710,
                    "1", "SARS-CoV-2_17_LEFT-SARS-CoV-2_17_RIGHT", 4923, 5331,
                    "1", "SARS-CoV-2_19_LEFT-SARS-CoV-2_19_RIGHT", 5561, 5957,
                    "1", "SARS-CoV-2_21_LEFT-SARS-CoV-2_21_RIGHT", 6184, 6582,
                    "1", "SARS-CoV-2_23_LEFT-SARS-CoV-2_23_RIGHT", 6747, 7148,
                    "1", "SARS-CoV-2_25_LEFT-SARS-CoV-2_25_RIGHT", 7381, 7770,
                    "1", "SARS-CoV-2_27_LEFT-SARS-CoV-2_27_RIGHT", 7997, 8395,
                    "1", "SARS-CoV-2_29_LEFT-SARS-CoV-2_29_RIGHT", 8596, 9013,
                    "1", "SARS-CoV-2_31_LEFT-SARS-CoV-2_31_RIGHT", 9168, 9564,
                    "1", "SARS-CoV-2_33_LEFT-SARS-CoV-2_33_RIGHT", 9782, 10176,
                    "1", "SARS-CoV-2_35_LEFT-SARS-CoV-2_35_RIGHT", 10393, 10810,
                    "1", "SARS-CoV-2_37_LEFT-SARS-CoV-2_37_RIGHT", 11000, 11414,
                    "1", "SARS-CoV-2_39_LEFT-SARS-CoV-2_39_RIGHT", 11624, 12033,
                    "1", "SARS-CoV-2_41_LEFT-SARS-CoV-2_41_RIGHT", 12234, 12643,
                    "1", "SARS-CoV-2_43_LEFT-SARS-CoV-2_43_RIGHT", 12831, 13240,
                    "1", "SARS-CoV-2_45_LEFT-SARS-CoV-2_45_RIGHT", 13463, 13859,
                    "1", "SARS-CoV-2_47_LEFT-SARS-CoV-2_47_RIGHT", 14045, 14457,
                    "1", "SARS-CoV-2_49_LEFT-SARS-CoV-2_49_RIGHT", 14647, 15050,
                    "1", "SARS-CoV-2_51_LEFT-SARS-CoV-2_51_RIGHT", 15214, 15619,
                    "1", "SARS-CoV-2_53_LEFT-SARS-CoV-2_53_RIGHT", 15855, 16260,
                    "1", "SARS-CoV-2_55_LEFT-SARS-CoV-2_55_RIGHT", 16386, 16796,
                    "1", "SARS-CoV-2_57_LEFT-SARS-CoV-2_57_RIGHT", 16986, 17405,
                    "1", "SARS-CoV-2_59_LEFT-SARS-CoV-2_59_RIGHT", 17615, 18022,
                    "1", "SARS-CoV-2_61_LEFT-SARS-CoV-2_61_RIGHT", 18244, 18652,
                    "1", "SARS-CoV-2_63_LEFT-SARS-CoV-2_63_RIGHT", 18869, 19277,
                    "1", "SARS-CoV-2_65_LEFT-SARS-CoV-2_65_RIGHT", 19485, 19901,
                    "1", "SARS-CoV-2_67_LEFT-SARS-CoV-2_67_RIGHT", 20090, 20497,
                    "1", "SARS-CoV-2_69_LEFT-SARS-CoV-2_69_RIGHT", 20677, 21080,
                    "1", "SARS-CoV-2_71_LEFT-SARS-CoV-2_71_RIGHT", 21294, 21700,
                    "1", "SARS-CoV-2_73_LEFT-SARS-CoV-2_73_RIGHT", 21865, 22274,
                    "1", "SARS-CoV-2_75_LEFT-SARS-CoV-2_75_RIGHT", 22402, 22805,
                    "1", "SARS-CoV-2_77_LEFT-SARS-CoV-2_77_RIGHT", 22944, 23351,
                    "1", "SARS-CoV-2_79_LEFT-SARS-CoV-2_79_RIGHT", 23553, 23955,
                    "1", "SARS-CoV-2_81_LEFT-SARS-CoV-2_81_RIGHT", 24171, 24567,
                    "1", "SARS-CoV-2_83_LEFT-SARS-CoV-2_83_RIGHT", 24750, 25150,
                    "1", "SARS-CoV-2_85_LEFT-SARS-CoV-2_85_RIGHT", 25331, 25740,
                    "1", "SARS-CoV-2_87_LEFT-SARS-CoV-2_87_RIGHT", 25951, 26360,
                    "1", "SARS-CoV-2_89_LEFT-SARS-CoV-2_89_RIGHT", 26564, 26979,
                    "1", "SARS-CoV-2_91_LEFT-SARS-CoV-2_91_RIGHT", 27152, 27560,
                    "1", "SARS-CoV-2_93_LEFT-SARS-CoV-2_93_RIGHT", 27700, 28104,
                    "1", "SARS-CoV-2_95_LEFT-SARS-CoV-2_95_RIGHT", 28190, 28598,
                    "1", "SARS-CoV-2_97_LEFT-SARS-CoV-2_97_RIGHT", 28827, 29227,
                    "1", "SARS-CoV-2_99_LEFT-SARS-CoV-2_99_RIGHT", 29452, 29854)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "SARS-CoV-2_23_LEFT-SARS-CoV-2_23_RIGHT_alt1", 6747, 7156,
                    "1", "SARS-CoV-2_27_LEFT-SARS-CoV-2_27_RIGHT_alt1", 7997, 8392,
                    "1", "SARS-CoV-2_79_LEFT-SARS-CoV-2_79_RIGHT_alt1", 23553, 23944,
                    "1", "SARS-CoV-2_89_LEFT_alt1-SARS-CoV-2_89_RIGHT_alt1", 26592, 26991)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p3 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "SARS-CoV-2_2_LEFT-SARS-CoV-2_2_RIGHT", 324, 727,
                    "2", "SARS-CoV-2_4_LEFT-SARS-CoV-2_4_RIGHT", 944, 1362,
                    "2", "SARS-CoV-2_6_LEFT-SARS-CoV-2_6_RIGHT", 1540, 1948,
                    "2", "SARS-CoV-2_8_LEFT-SARS-CoV-2_8_RIGHT", 2154, 2571,
                    "2", "SARS-CoV-2_10_LEFT-SARS-CoV-2_10_RIGHT", 2826, 3210,
                    "2", "SARS-CoV-2_12_LEFT-SARS-CoV-2_12_RIGHT", 3390, 3794,
                    "2", "SARS-CoV-2_14_LEFT-SARS-CoV-2_14_RIGHT", 3992, 4409,
                    "2", "SARS-CoV-2_16_LEFT-SARS-CoV-2_16_RIGHT", 4620, 5017,
                    "2", "SARS-CoV-2_18_LEFT-SARS-CoV-2_18_RIGHT", 5230, 5643,
                    "2", "SARS-CoV-2_20_LEFT-SARS-CoV-2_20_RIGHT", 5867, 6272,
                    "2", "SARS-CoV-2_22_LEFT-SARS-CoV-2_22_RIGHT", 6478, 6885,
                    "2", "SARS-CoV-2_24_LEFT-SARS-CoV-2_24_RIGHT", 7057, 7467,
                    "2", "SARS-CoV-2_26_LEFT-SARS-CoV-2_26_RIGHT", 7672, 8092,
                    "2", "SARS-CoV-2_28_LEFT-SARS-CoV-2_28_RIGHT", 8304, 8714,
                    "2", "SARS-CoV-2_30_LEFT-SARS-CoV-2_30_RIGHT", 8919, 9329,
                    "2", "SARS-CoV-2_32_LEFT-SARS-CoV-2_32_RIGHT", 9470, 9866,
                    "2", "SARS-CoV-2_34_LEFT-SARS-CoV-2_34_RIGHT", 10076, 10491,
                    "2", "SARS-CoV-2_36_LEFT-SARS-CoV-2_36_RIGHT", 10713, 11116,
                    "2", "SARS-CoV-2_38_LEFT-SARS-CoV-2_38_RIGHT", 11305, 11720,
                    "2", "SARS-CoV-2_40_LEFT-SARS-CoV-2_40_RIGHT", 11937, 12339,
                    "2", "SARS-CoV-2_42_LEFT-SARS-CoV-2_42_RIGHT", 12519, 12920,
                    "2", "SARS-CoV-2_44_LEFT-SARS-CoV-2_44_RIGHT", 13124, 13528,
                    "2", "SARS-CoV-2_46_LEFT-SARS-CoV-2_46_RIGHT", 13752, 14144,
                    "2", "SARS-CoV-2_48_LEFT-SARS-CoV-2_48_RIGHT", 14338, 14743,
                    "2", "SARS-CoV-2_50_LEFT-SARS-CoV-2_50_RIGHT", 14953, 15358,
                    "2", "SARS-CoV-2_52_LEFT-SARS-CoV-2_52_RIGHT", 15535, 15941,
                    "2", "SARS-CoV-2_54_LEFT-SARS-CoV-2_54_RIGHT", 16112, 16508,
                    "2", "SARS-CoV-2_56_LEFT-SARS-CoV-2_56_RIGHT", 16692, 17105,
                    "2", "SARS-CoV-2_58_LEFT-SARS-CoV-2_58_RIGHT", 17323, 17711,
                    "2", "SARS-CoV-2_60_LEFT-SARS-CoV-2_60_RIGHT", 17911, 18328,
                    "2", "SARS-CoV-2_62_LEFT-SARS-CoV-2_62_RIGHT", 18550, 18961,
                    "2", "SARS-CoV-2_64_LEFT-SARS-CoV-2_64_RIGHT", 19183, 19586,
                    "2", "SARS-CoV-2_66_LEFT-SARS-CoV-2_66_RIGHT", 19810, 20216,
                    "2", "SARS-CoV-2_68_LEFT-SARS-CoV-2_68_RIGHT", 20377, 20792,
                    "2", "SARS-CoV-2_70_LEFT-SARS-CoV-2_70_RIGHT", 20988, 21387,
                    "2", "SARS-CoV-2_72_LEFT-SARS-CoV-2_72_RIGHT", 21532, 21933,
                    "2", "SARS-CoV-2_74_LEFT-SARS-CoV-2_74_RIGHT", 22091, 22503,
                    "2", "SARS-CoV-2_76_LEFT-SARS-CoV-2_76_RIGHT", 22648, 23057,
                    "2", "SARS-CoV-2_78_LEFT-SARS-CoV-2_78_RIGHT", 23219, 23635,
                    "2", "SARS-CoV-2_80_LEFT-SARS-CoV-2_80_RIGHT", 23853, 24258,
                    "2", "SARS-CoV-2_82_LEFT-SARS-CoV-2_82_RIGHT", 24426, 24836,
                    "2", "SARS-CoV-2_84_LEFT-SARS-CoV-2_84_RIGHT", 25051, 25461,
                    "2", "SARS-CoV-2_86_LEFT-SARS-CoV-2_86_RIGHT", 25645, 26050,
                    "2", "SARS-CoV-2_88_LEFT-SARS-CoV-2_88_RIGHT", 26255, 26661,
                    "2", "SARS-CoV-2_90_LEFT-SARS-CoV-2_90_RIGHT", 26873, 27283,
                    "2", "SARS-CoV-2_92_LEFT-SARS-CoV-2_92_RIGHT", 27447, 27855,
                    "2", "SARS-CoV-2_94_LEFT-SARS-CoV-2_94_RIGHT", 27996, 28416,
                    "2", "SARS-CoV-2_96_LEFT-SARS-CoV-2_96_RIGHT", 28512, 28914,
                    "2", "SARS-CoV-2_98_LEFT-SARS-CoV-2_98_RIGHT", 29136, 29534)
  map1plot3 <- map1p3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map1p4 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "SARS-CoV-2_10_LEFT_alt1-SARS-CoV-2_10_RIGHT_alt1", 2780, 3177,
                    "2", "SARS-CoV-2_76_LEFT_alt1-SARS-CoV-2_76_RIGHT_alt1", 22742, 23141,
                    "2", "SARS-CoV-2_88_LEFT_alt1-SARS-CoV-2_88_RIGHT", 26242, 26661,
                    "2", "SARS-CoV-2_90_LEFT-SARS-CoV-2_90_RIGHT_alt1", 26873, 27251)
  map1plot4 <- map1p4 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / map1plot3 / map1plot4 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 9, heights = c(3, .1, .1, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ARTIC/V5.3.2") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "SARS-CoV-2_400_1_LEFT_1-SARS-CoV-2_400_1_RIGHT_1", 47, 447,
                    "1", "SARS-CoV-2_400_3_LEFT_1-SARS-CoV-2_400_3_RIGHT_0", 638, 1047,
                    "1", "SARS-CoV-2_400_5_LEFT_0-SARS-CoV-2_400_5_RIGHT_0", 1292, 1692,
                    "1", "SARS-CoV-2_400_7_LEFT_2-SARS-CoV-2_400_7_RIGHT_2", 1882, 2284,
                    "1", "SARS-CoV-2_400_9_LEFT_0-SARS-CoV-2_400_9_RIGHT_0", 2533, 2933,
                    "1", "SARS-CoV-2_400_11_LEFT_0-SARS-CoV-2_400_11_RIGHT_0", 3184, 3584,
                    "1", "SARS-CoV-2_400_13_LEFT_0-SARS-CoV-2_400_13_RIGHT_0", 3791, 4180,
                    "1", "SARS-CoV-2_400_15_LEFT_0-SARS-CoV-2_400_15_RIGHT_0", 4403, 4803,
                    "1", "SARS-CoV-2_400_17_LEFT_0-SARS-CoV-2_400_17_RIGHT_0", 5036, 5429,
                    "1", "SARS-CoV-2_400_19_LEFT_0-SARS-CoV-2_400_19_RIGHT_0", 5671, 6062,
                    "1", "SARS-CoV-2_400_21_LEFT_0-SARS-CoV-2_400_21_RIGHT_0", 6204, 6595,
                    "1", "SARS-CoV-2_400_23_LEFT_0-SARS-CoV-2_400_23_RIGHT_0", 6823, 7229,
                    "1", "SARS-CoV-2_400_25_LEFT_0-SARS-CoV-2_400_25_RIGHT_0", 7456, 7850,
                    "1", "SARS-CoV-2_400_27_LEFT_0-SARS-CoV-2_400_27_RIGHT_0", 8085, 8498,
                    "1", "SARS-CoV-2_400_29_LEFT_0-SARS-CoV-2_400_29_RIGHT_0", 8732, 9129,
                    "1", "SARS-CoV-2_400_31_LEFT_1-SARS-CoV-2_400_31_RIGHT_0", 9299, 9706,
                    "1", "SARS-CoV-2_400_33_LEFT_0-SARS-CoV-2_400_33_RIGHT_0", 9896, 10295,
                    "1", "SARS-CoV-2_400_35_LEFT_0-SARS-CoV-2_400_35_RIGHT_0", 10527, 10927,
                    "1", "SARS-CoV-2_400_37_LEFT_0-SARS-CoV-2_400_37_RIGHT_0", 11152, 11536,
                    "1", "SARS-CoV-2_400_39_LEFT_0-SARS-CoV-2_400_39_RIGHT_0", 11785, 12185,
                    "1", "SARS-CoV-2_400_41_LEFT_0-SARS-CoV-2_400_41_RIGHT_0", 12419, 12819,
                    "1", "SARS-CoV-2_400_43_LEFT_0-SARS-CoV-2_400_43_RIGHT_0", 13075, 13480,
                    "1", "SARS-CoV-2_400_45_LEFT_0-SARS-CoV-2_400_45_RIGHT_0", 13738, 14144,
                    "1", "SARS-CoV-2_400_47_LEFT_0-SARS-CoV-2_400_47_RIGHT_0", 14375, 14775,
                    "1", "SARS-CoV-2_400_49_LEFT_0-SARS-CoV-2_400_49_RIGHT_0", 15016, 15416,
                    "1", "SARS-CoV-2_400_51_LEFT_0-SARS-CoV-2_400_51_RIGHT_0", 15659, 16059,
                    "1", "SARS-CoV-2_400_53_LEFT_0-SARS-CoV-2_400_53_RIGHT_0", 16285, 16679,
                    "1", "SARS-CoV-2_400_55_LEFT_1-SARS-CoV-2_400_55_RIGHT_1", 16962, 17362,
                    "1", "SARS-CoV-2_400_57_LEFT_0-SARS-CoV-2_400_57_RIGHT_0", 17478, 17886,
                    "1", "SARS-CoV-2_400_59_LEFT_0-SARS-CoV-2_400_59_RIGHT_0", 18121, 18527,
                    "1", "SARS-CoV-2_400_61_LEFT_0-SARS-CoV-2_400_61_RIGHT_0", 18789, 19195,
                    "1", "SARS-CoV-2_400_63_LEFT_0-SARS-CoV-2_400_63_RIGHT_0", 19415, 19796,
                    "1", "SARS-CoV-2_400_65_LEFT_0-SARS-CoV-2_400_65_RIGHT_0", 20028, 20441,
                    "1", "SARS-CoV-2_400_67_LEFT_1-SARS-CoV-2_400_67_RIGHT_1", 20650, 21051,
                    "1", "SARS-CoV-2_400_69_LEFT_0-SARS-CoV-2_400_69_RIGHT_0", 21322, 21722,
                    "1", "SARS-CoV-2_400_71_LEFT_0-SARS-CoV-2_400_71_RIGHT_0", 21866, 22266,
                    "1", "SARS-CoV-2_400_73_LEFT_0-SARS-CoV-2_400_73_RIGHT_0", 22466, 22866,
                    "1", "SARS-CoV-2_400_75_LEFT_1-SARS-CoV-2_400_75_RIGHT_1", 23078, 23478,
                    "1", "SARS-CoV-2_400_77_LEFT_0-SARS-CoV-2_400_77_RIGHT_0", 23563, 23944,
                    "1", "SARS-CoV-2_400_79_LEFT_0-SARS-CoV-2_400_79_RIGHT_0", 24160, 24560,
                    "1", "SARS-CoV-2_400_81_LEFT_0-SARS-CoV-2_400_81_RIGHT_0", 24751, 25151,
                    "1", "SARS-CoV-2_400_83_LEFT_0-SARS-CoV-2_400_83_RIGHT_0", 25372, 25777,
                    "1", "SARS-CoV-2_400_85_LEFT_0-SARS-CoV-2_400_85_RIGHT_0", 26011, 26411,
                    "1", "SARS-CoV-2_400_87_LEFT_1-SARS-CoV-2_400_87_RIGHT_1", 26593, 27009,
                    "1", "SARS-CoV-2_400_89_LEFT_2-SARS-CoV-2_400_89_RIGHT_0", 27200, 27603,
                    "1", "SARS-CoV-2_400_91_LEFT_0-SARS-CoV-2_400_91_RIGHT_0", 27832, 28237,
                    "1", "SARS-CoV-2_400_93_LEFT_0-SARS-CoV-2_400_93_RIGHT_0", 28473, 28873,
                    "1", "SARS-CoV-2_400_95_LEFT_0-SARS-CoV-2_400_95_RIGHT_0", 29159, 29559)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "SARS-CoV-2_400_2_LEFT_0-SARS-CoV-2_400_2_RIGHT_0", 344, 732,
                    "2", "SARS-CoV-2_400_4_LEFT_0-SARS-CoV-2_400_4_RIGHT_0", 970, 1370,
                    "2", "SARS-CoV-2_400_6_LEFT_1-SARS-CoV-2_400_6_RIGHT_1", 1574, 1972,
                    "2", "SARS-CoV-2_400_8_LEFT_0-SARS-CoV-2_400_8_RIGHT_0", 2229, 2629,
                    "2", "SARS-CoV-2_400_10_LEFT_0-SARS-CoV-2_400_10_RIGHT_0", 2854, 3254,
                    "2", "SARS-CoV-2_400_12_LEFT_0-SARS-CoV-2_400_12_RIGHT_0", 3510, 3913,
                    "2", "SARS-CoV-2_400_14_LEFT_0-SARS-CoV-2_400_14_RIGHT_0", 4079, 4488,
                    "2", "SARS-CoV-2_400_16_LEFT_0-SARS-CoV-2_400_16_RIGHT_0", 4723, 5119,
                    "2", "SARS-CoV-2_400_18_LEFT_0-SARS-CoV-2_400_18_RIGHT_0", 5344, 5744,
                    "2", "SARS-CoV-2_400_20_LEFT_0-SARS-CoV-2_400_20_RIGHT_0", 5891, 6288,
                    "2", "SARS-CoV-2_400_22_LEFT_0-SARS-CoV-2_400_22_RIGHT_0", 6515, 6915,
                    "2", "SARS-CoV-2_400_24_LEFT_0-SARS-CoV-2_400_24_RIGHT_0", 7145, 7545,
                    "2", "SARS-CoV-2_400_26_LEFT_0-SARS-CoV-2_400_26_RIGHT_0", 7768, 8169,
                    "2", "SARS-CoV-2_400_28_LEFT_0-SARS-CoV-2_400_28_RIGHT_0", 8406, 8806,
                    "2", "SARS-CoV-2_400_30_LEFT_0-SARS-CoV-2_400_30_RIGHT_0", 9023, 9423,
                    "2", "SARS-CoV-2_400_32_LEFT_0-SARS-CoV-2_400_32_RIGHT_0", 9571, 9971,
                    "2", "SARS-CoV-2_400_34_LEFT_0-SARS-CoV-2_400_34_RIGHT_0", 10215, 10615,
                    "2", "SARS-CoV-2_400_36_LEFT_0-SARS-CoV-2_400_36_RIGHT_0", 10832, 11232,
                    "2", "SARS-CoV-2_400_38_LEFT_0-SARS-CoV-2_400_38_RIGHT_0", 11463, 11863,
                    "2", "SARS-CoV-2_400_40_LEFT_0-SARS-CoV-2_400_40_RIGHT_0", 12112, 12510,
                    "2", "SARS-CoV-2_400_42_LEFT_0-SARS-CoV-2_400_42_RIGHT_0", 12752, 13146,
                    "2", "SARS-CoV-2_400_44_LEFT_0-SARS-CoV-2_400_44_RIGHT_0", 13415, 13815,
                    "2", "SARS-CoV-2_400_46_LEFT_0-SARS-CoV-2_400_46_RIGHT_0", 14073, 14457,
                    "2", "SARS-CoV-2_400_48_LEFT_0-SARS-CoV-2_400_48_RIGHT_0", 14700, 15095,
                    "2", "SARS-CoV-2_400_50_LEFT_0-SARS-CoV-2_400_50_RIGHT_0", 15342, 15742,
                    "2", "SARS-CoV-2_400_52_LEFT_2-SARS-CoV-2_400_52_RIGHT_2", 15992, 16409,
                    "2", "SARS-CoV-2_400_54_LEFT_1-SARS-CoV-2_400_54_RIGHT_1", 16624, 17033,
                    "2", "SARS-CoV-2_400_56_LEFT_0-SARS-CoV-2_400_56_RIGHT_0", 17182, 17582,
                    "2", "SARS-CoV-2_400_58_LEFT_0-SARS-CoV-2_400_58_RIGHT_0", 17813, 18212,
                    "2", "SARS-CoV-2_400_60_LEFT_0-SARS-CoV-2_400_60_RIGHT_0", 18460, 18860,
                    "2", "SARS-CoV-2_400_62_LEFT_2-SARS-CoV-2_400_62_RIGHT_0", 19087, 19495,
                    "2", "SARS-CoV-2_400_64_LEFT_0-SARS-CoV-2_400_64_RIGHT_0", 19721, 20121,
                    "2", "SARS-CoV-2_400_66_LEFT_0-SARS-CoV-2_400_66_RIGHT_0", 20358, 20758,
                    "2", "SARS-CoV-2_400_68_LEFT_0-SARS-CoV-2_400_68_RIGHT_0", 20991, 21402,
                    "2", "SARS-CoV-2_400_70_LEFT_0-SARS-CoV-2_400_70_RIGHT_0", 21579, 21960,
                    "2", "SARS-CoV-2_400_72_LEFT_0-SARS-CoV-2_400_72_RIGHT_0", 22156, 22547,
                    "2", "SARS-CoV-2_400_74_LEFT_0-SARS-CoV-2_400_74_RIGHT_0", 22742, 23140,
                    "2", "SARS-CoV-2_400_76_LEFT_0-SARS-CoV-2_400_76_RIGHT_0", 23229, 23631,
                    "2", "SARS-CoV-2_400_78_LEFT_0-SARS-CoV-2_400_78_RIGHT_0", 23823, 24231,
                    "2", "SARS-CoV-2_400_80_LEFT_0-SARS-CoV-2_400_80_RIGHT_0", 24442, 24839,
                    "2", "SARS-CoV-2_400_82_LEFT_0-SARS-CoV-2_400_82_RIGHT_0", 25053, 25452,
                    "2", "SARS-CoV-2_400_84_LEFT_2-SARS-CoV-2_400_84_RIGHT_2", 25653, 26072,
                    "2", "SARS-CoV-2_400_86_LEFT_0-SARS-CoV-2_400_86_RIGHT_0", 26339, 26756,
                    "2", "SARS-CoV-2_400_88_LEFT_2-SARS-CoV-2_400_88_RIGHT_2", 26958, 27376,
                    "2", "SARS-CoV-2_400_90_LEFT_0-SARS-CoV-2_400_90_RIGHT_0", 27530, 27950,
                    "2", "SARS-CoV-2_400_92_LEFT_0-SARS-CoV-2_400_92_RIGHT_0", 28135, 28539,
                    "2", "SARS-CoV-2_400_94_LEFT_0-SARS-CoV-2_400_94_RIGHT_0", 28808, 29224,
                    "2", "SARS-CoV-2_400_96_LEFT_1-SARS-CoV-2_400_96_RIGHT_0", 29462, 29873)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = .5) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "FIOCRUZ-IOC/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "hCoV_F1_2kb-hCoV_R1_2kb", 31, 2079,
                    "1", "hCoV_F3_2kb-hCoV_R3_2kb", 3581, 5548,
                    "1", "hCoV_F5_2kb-hCoV_R5_2kb", 7093, 9123,
                    "1", "hCoV_F7_2kb-hCoV_R7_2kb", 10677, 12679,
                    "1", "hCoV_F9_2kb-hCoV_R9_2kb", 17177, 15978,
                    "1", "hCoV_F11_2kb-hCoV_R11_2kb", 17572, 19485,
                    "1", "hCoV_F13_2kb-hCoV_R13_2kb", 21076, 22996,
                    "1", "hCoV_F15_2kb-hCoV_R15_2kb", 24650, 26542,
                    "1", "hCoV_F17_2kb-hCoV_R17_2kb", 27915, 29790)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "hCoV_F2_2kb-hCoV_R2_2kb", 1926, 3737,
                    "2", "hCoV_F4_2kb-hCoV_R4_2kb", 5395, 7255,
                    "2", "hCoV_F6_2kb-hCoV_R6_2kb", 8799, 10837,
                    "2", "hCoV_F8_2kb-hCoV_R8_2kb", 12520, 14328,
                    "2", "hCoV_F10_2kb-hCoV_R10_2kb", 15828, 17754,
                    "2", "hCoV_F12_2kb-hCoV_R12_2kb", 19311, 21241,
                    "2", "hCoV_F14_2kb-hCoV_R14_2kb", 22851, 24812,
                    "2", "hCoV_F16_2kb-hCoV_R16_2kb", 26387, 28351)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "FIOCRUZ-IOC/V2") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 5000, 10000, 15000, 20000, 25000, 29903),
                       expand = expansion(0, 0), limits = c(0, 30500)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "hCoV_F1_2kb-hCoV_R1_2kb", 31, 2079,
                    "1", "hCoV_F3_2kb-hCoV_R3_2kb", 3581, 5548,
                    "1", "hCoV_F5_2kb-hCoV_R5_2kb", 7093, 9123,
                    "1", "hCoV_F7_2kb-hCoV_R7_2kb", 10677, 12679,
                    "1", "hCoV_F9_2kb-hCoV_R9_2kb", 14177, 15978,
                    "1", "hCoV_F11_2kb-hCoV_R11_2kb", 17572, 19485,
                    "1", "hCoV_F13_2kb-hCoV_R13_2kb", 21076, 22996,
                    "1", "hCoV_F15_2kb-hCoV_R15_2kb", 24650, 26542,
                    "1", "hCoV_F17_2kb-hCoV_R17_2kb", 27915, 29790)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "hCoV_F2_2kb-hCoV_R2_2kb", 1926, 3737,
                    "2", "hCoV_F4_2kb-hCoV_R4_2kb", 5395, 7255,
                    "2", "hCoV_F6_2kb-hCoV_R6_2kb", 8799, 10837,
                    "2", "hCoV_F8_2kb-hCoV_R8_2kb", 12520, 14328,
                    "2", "hCoV_F10_2kb-hCoV_R10_2kb", 15828, 17754,
                    "2", "hCoV_F12_2kb-hCoV_R12_2kb", 19311, 21241,
                    "2", "hCoV_F14_2kb-hCoV_R14_2kb", 22851, 24812,
                    "2", "hCoV_F16_2kb-hCoV_R16_2kb", 26387, 28351)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map1p3 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "hCoV_F44_2kb-hCoV_R8_2kb", 12160, 14328,
                    "2", "hCoV_F55_2kb-hCoV_R10_2kb", 15375, 17754,)
  map1plot3 <- map1p3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "5'UTR", 1, 265,
                         "ORF1a", 266, 13468,
                         "ORF1b", 13468, 21555,
                         "S", 21563, 25384,
                         "E", 26245, 26472,
                         "M", 26523, 27191,
                         "N", 28274, 29533,
                         "3'UTR", 29675, 29903)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp1", 266, 805,
                         "nsp2", 806, 2719,
                         "nsp3", 2720, 8554,
                         "nsp4", 8555, 10054,
                         "nsp5", 10055, 10972,
                         "nsp6", 10973, 11842,
                         "nsp8", 12092, 12685,
                         "nsp10", 13025, 13441,
                         "nsp12", 13442, 16236,
                         "nsp13", 16237, 18039,
                         "nsp14", 18040, 19620,
                         "nsp15", 19621, 20658,
                         "nsp16", 20659, 21552,
                         "3ab", 25393, 26220,
                         "6", 27202, 27387,
                         "7a", 27394, 27759,
                         "7b", 27756, 27887,
                         "8", 27894, 28259,
                         "9b", 28284, 28577,
                         "10", 29558, 28674)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/1798174254
                                                    # https://doi.org/10.1038/s41586-020-2286-9
                         "nsp7", 11843, 12091,
                         "nsp9", 12686, 13024,
                         "nsp11", 13442, 13480)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 30500)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".sars2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / map1plot3 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 8, heights = c(3, .1, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ZikaAsian/V2") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10807),
                       expand = expansion(0, 0), limits = c(0, 10810)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "ZIKA_400_1_L-ZIKA_400_1_R", 28, 504,
                    "1", "ZIKA_400_3_L-ZIKA_400_3_R", 658, 1076,
                    "1", "ZIKA_400_5_L-ZIKA_400_5_R", 1256, 1723,
                    "1", "ZIKA_400_7_L-ZIKA_400_7_R", 1875, 2286,
                    "1", "ZIKA_400_9_L-ZIKA_400_9_R", 2444, 2924,
                    "1", "ZIKA_400_11_L-ZIKA_400_11_R", 3030, 3503,
                    "1", "ZIKA_400_13_L-ZIKA_400_13_R", 3657, 4099,
                    "1", "ZIKA_400_15_L-ZIKA_400_15_R", 4231, 4681,
                    "1", "ZIKA_400_17_L-ZIKA_400_17_R", 4851, 5326,
                    "1", "ZIKA_400_19_L-ZIKA_400_19_R", 5461, 5924,
                    "1", "ZIKA_400_21_L-ZIKA_400_21_R", 6046, 6480,
                    "1", "ZIKA_400_23_L-ZIKA_400_23_R", 6673, 7095,
                    "1", "ZIKA_400_25_L-ZIKA_400_25_R", 7230, 7707,
                    "1", "ZIKA_400_27_L-ZIKA_400_27_R", 7841, 8318,
                    "1", "ZIKA_400_29_L-ZIKA_400_29_R", 8430, 8895,
                    "1", "ZIKA_400_31_L-ZIKA_400_31_R", 9052, 9495,
                    "1", "ZIKA_400_33_L-ZIKA_400_33_R", 9637, 10126,
                    "1", "ZIKA_400_35_L-ZIKA_400_35_R", 10241, 10669,)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10810)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "ZIKA_400_2_L-ZIKA_400_2_R", 359, 818,
                    "2", "ZIKA_400_4_L-ZIKA_400_4_R", 978, 1415,
                    "2", "ZIKA_400_6_L-ZIKA_400_6_R", 1539, 1987,
                    "2", "ZIKA_400_8_L-ZIKA_400_8_R", 2151, 2629,
                    "2", "ZIKA_400_10_L-ZIKA_400_10_R", 2747, 3176,
                    "2", "ZIKA_400_12_L-ZIKA_400_12_R", 3334, 3785,
                    "2", "ZIKA_400_14_L-ZIKA_400_14_R", 3931, 4410,
                    "2", "ZIKA_400_16_L-ZIKA_400_16_R", 4542, 4987,
                    "2", "ZIKA_400_18_L-ZIKA_400_18_R", 5136, 5609,
                    "2", "ZIKA_400_20_L-ZIKA_400_20_R", 5746, 6200,
                    "2", "ZIKA_400_22_L-ZIKA_400_22_R", 6353, 6778,
                    "2", "ZIKA_400_24_L-ZIKA_400_24_R", 6960, 7416,
                    "2", "ZIKA_400_26_L-ZIKA_400_26_R", 7554, 7994,
                    "2", "ZIKA_400_28_L-ZIKA_400_28_R", 8167, 8604,
                    "2", "ZIKA_400_30_L-ZIKA_400_30_R", 8740, 9206,
                    "2", "ZIKA_400_32_L-ZIKA_400_32_R", 9336, 9777,
                    "2", "ZIKA_400_34_L-ZIKA_400_34_R", 9969, 10413)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10810)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", #  https://www.ncbi.nlm.nih.gov/nuccore/KJ776791.2
                                                    # https://doi.org/10.1371/journal.ppat.1006528
                         "5'UTR", 1, 107,
                         "pr", 474, 752,
                         "M", 753, 977,
                         "E", 978, 2489,
                         "NS1", 2490, 3545,
                         "NS2A", 3546, 4223,
                         "NS3", 4614, 6464,
                         "NS4A", 6465, 6845,
                         "NS4B", 6915, 7667,
                         "NS5", 7668, 10379,
                         "3'UTR", 10380, 10807)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10810)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", #  https://www.ncbi.nlm.nih.gov/nuccore/KJ776791.2
                                                    # https://doi.org/10.1371/journal.ppat.1006528
                         "C", 108, 473,
                         "NS2B", 4224, 4613,
                         "2K", 6846, 6914)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10810)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".zikv-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENGUESEQ1/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10735),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV1_1_LEFT-DENV1_1_RIGHT", 12, 477,
                    "1", "DENV1_3_LEFT-DENV1_3_RIGHT", 691, 1093,
                    "1", "DENV1_5_LEFT-DENV1_5_RIGHT", 1303, 1720,
                    "1", "DENV1_7_LEFT-DENV1_7_RIGHT", 1912, 2329,
                    "1", "DENV1_9_LEFT-DENV1_9_RIGHT", 2487, 2883,
                    "1", "DENV1_11_LEFT-DENV1_11_RIGHT", 3101, 3489,
                    "1", "DENV1_13_LEFT-DENV1_13_RIGHT", 3641, 4039,
                    "1", "DENV1_15_LEFT-DENV1_15_RIGHT", 4114, 4524,
                    "1", "DENV1_17_LEFT-DENV1_17_RIGHT", 4694, 5111,
                    "1", "DENV1_19_LEFT-DENV1_19_RIGHT", 5332, 5742,
                    "1", "DENV1_21_LEFT-DENV1_21_RIGHT", 5974, 6376,
                    "1", "DENV1_23_LEFT-DENV1_23_RIGHT", 6475, 6868,
                    "1", "DENV1_25_LEFT-DENV1_25_RIGHT", 7089, 7482,
                    "1", "DENV1_27_LEFT-DENV1_27_RIGHT", 7710, 8097,
                    "1", "DENV1_29_LEFT-DENV1_29_RIGHT", 8299, 8709,
                    "1", "DENV1_31_LEFT-DENV1_31_RIGHT", 8917, 9326,
                    "1", "DENV1_33_LEFT-DENV1_33_RIGHT", 9533, 9925,
                    "1", "DENV1_35_LEFT-DENV1_35_RIGHT", 10143, 10531)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV1_2_LEFT-DENV1_2_RIGHT", 398, 786,
                    "2", "DENV1_4_LEFT-DENV1_4_RIGHT", 1015, 1410,
                    "2", "DENV1_6_LEFT-DENV1_6_RIGHT", 1640, 2033,
                    "2", "DENV1_8_LEFT-DENV1_8_RIGHT", 2241, 2640,
                    "2", "DENV1_10_LEFT-DENV1_10_RIGHT", 2763, 3182,
                    "2", "DENV1_12_LEFT-DENV1_12_RIGHT", 3328, 3735,
                    "2", "DENV1_14_LEFT-DENV1_14_RIGHT", 3856, 4248,
                    "2", "DENV1_16_LEFT-DENV1_16_RIGHT", 4381, 4800,
                    "2", "DENV1_18_LEFT-DENV1_18_RIGHT", 5014, 5418,
                    "2", "DENV1_20_LEFT-DENV1_20_RIGHT", 5659, 6053,
                    "2", "DENV1_22_LEFT-DENV1_22_RIGHT", 6194, 6582,
                    "2", "DENV1_24_LEFT-DENV1_24_RIGHT", 6774, 7168,
                    "2", "DENV1_26_LEFT-DENV1_26_RIGHT", 7389, 7791,
                    "2", "DENV1_28_LEFT-DENV1_28_RIGHT", 8004, 8408,
                    "2", "DENV1_30_LEFT-DENV1_30_RIGHT", 8606, 9012,
                    "2", "DENV1_32_LEFT-DENV1_32_RIGHT", 9212, 9614,
                    "2", "DENV1_34_LEFT-DENV1_34_RIGHT", 9841, 10241)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1
                         "5'UTR", 1, 94,
                         "prM", 437, 934,
                         "E", 935, 2419,
                         "NS1", 2420, 3475,
                         "NS2A", 3476, 4129,
                         "NS3", 4520, 6376,
                         "NS4A", 6377, 6757,
                         "NS4B", 6827, 7573,
                         "NS5", 7574, 10270,
                         "3'UTR", 10274, 10735)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1
                         "C", 95, 436,
                         "NS2B", 4130, 4519,
                         "2K", 6758, 6826)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv1-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENGUESEQ2/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV2_1_LEFT2-DENV2_1_RIGHT", 13, 481,
                    "1", "DENV2_3_LEFT-DENV2_3_RIGHT", 648, 1046,
                    "1", "DENV2_5_LEFT-DENV2_5_RIGHT", 1192, 1589,
                    "1", "DENV2_7_LEFT-DENV2_7_RIGHT", 1757, 2150,
                    "1", "DENV2_9_LEFT-DENV2_9_RIGHT", 2324, 2727,
                    "1", "DENV2_11_LEFT-DENV2_11_RIGHT", 2850, 3250,
                    "1", "DENV2_13_LEFT-DENV2_13_RIGHT", 3387, 3792,
                    "1", "DENV2_15_LEFT-DENV2_15_RIGHT", 3989, 4393,
                    "1", "DENV2_17_LEFT-DENV2_17_RIGHT", 4467, 4859,
                    "1", "DENV2_19_LEFT-DENV2_19_RIGHT", 5007, 5426,
                    "1", "DENV2_21_LEFT-DENV2_21_RIGHT", 5568, 5980,
                    "1", "DENV2_23_LEFT-DENV2_23_RIGHT", 6119, 6522,
                    "1", "DENV2_25_LEFT-DENV2_25_RIGHT", 6666, 7055,
                    "1", "DENV2_27_LEFT-DENV2_27_RIGHT", 7240, 7643,
                    "1", "DENV2_29_LEFT-DENV2_29_RIGHT", 7822, 8231,
                    "1", "DENV2_31_LEFT-DENV2_31_RIGHT", 8424, 8821,
                    "1", "DENV2_33_LEFT-DENV2_33_RIGHT", 8967, 9368,
                    "1", "DENV2_35_LEFT-DENV2_35_RIGHT", 9532, 9950,
                    "1", "DENV2_37_LEFT-DENV2_37_RIGHT", 10176, 10567)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV2_2_LEFT-DENV2_2_RIGHT", 368, 777,
                    "2", "DENV2_4_LEFT-DENV2_4_RIGHT", 919, 1316,
                    "2", "DENV2_6_LEFT-DENV2_6_RIGHT", 1503, 1898,
                    "2", "DENV2_8_LEFT-DENV2_8_RIGHT", 2036, 2443,
                    "2", "DENV2_10_LEFT-DENV2_10_RIGHT", 2620, 3020,
                    "2", "DENV2_12_LEFT-DENV2_12_RIGHT", 3165, 3571,
                    "2", "DENV2_14_LEFT-DENV2_14_RIGHT", 3696, 4097,
                    "2", "DENV2_16_LEFT-DENV2_16_RIGHT", 4119, 4522,
                    "2", "DENV2_18_LEFT-DENV2_18_RIGHT", 4768, 5158,
                    "2", "DENV2_20_LEFT-DENV2_20_RIGHT", 5316, 5736,
                    "2", "DENV2_22_LEFT-DENV2_22_RIGHT", 5882, 6286,
                    "2", "DENV2_24_LEFT-DENV2_24_RIGHT", 6352, 6758,
                    "2", "DENV2_26_LEFT-DENV2_26_RIGHT", 6914, 7325,
                    "2", "DENV2_28_LEFT-DENV2_28_RIGHT", 7516, 7924,
                    "2", "DENV2_30_LEFT-DENV2_30_RIGHT", 8092, 8507,
                    "2", "DENV2_32_LEFT-DENV2_32_RIGHT", 8655, 9050,
                    "2", "DENV2_34_LEFT-DENV2_34_RIGHT", 9228, 9630,
                    "2", "DENV2_36_LEFT-DENV2_36_RIGHT", 9862, 10259)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001474.2
                         "5'UTR", 1, 96,
                         "prM", 439, 936,
                         "E", 937, 2421,
                         "NS1", 2422, 3477,
                         "NS2A", 3478, 4131,
                         "NS3", 4522, 6375,
                         "NS4A", 6376, 6756,
                         "NS4B", 6826, 7569,
                         "NS5", 7570, 10269,
                         "3'UTR", 10273, 10723)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001474.2
                         "C", 97, 438,
                         "NS2B", 4132, 4521,
                         "2K", 6757, 6825)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENGUESEQ3/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV3_1_LEFT-DENV3_1_RIGHT", 13, 416,
                    "1", "DENV3_3_LEFT-DENV3_3_RIGHT", 616, 1015,
                    "1", "DENV3_5_LEFT-DENV3_5_RIGHT", 1236, 1641,
                    "1", "DENV3_7_LEFT-DENV3_7_RIGHT", 1835, 2241,
                    "1", "DENV3_9_LEFT-DENV3_9_RIGHT", 2414, 2828,
                    "1", "DENV3_11_LEFT-DENV3_11_RIGHT", 3048, 3457,
                    "1", "DENV3_13_LEFT-DENV3_13_RIGHT", 3641, 4044,
                    "1", "DENV3_15_LEFT-DENV3_15_RIGHT", 4246, 4645,
                    "1", "DENV3_17_LEFT-DENV3_17_RIGHT", 4873, 5275,
                    "1", "DENV3_19_LEFT-DENV3_19_RIGHT", 5482, 5873,
                    "1", "DENV3_21_LEFT-DENV3_21_RIGHT", 6088, 6504,
                    "1", "DENV3_23_LEFT-DENV3_23_RIGHT", 6746, 7133,
                    "1", "DENV3_25_LEFT-DENV3_25_RIGHT", 7357, 7764,
                    "1", "DENV3_27_LEFT-DENV3_27_RIGHT", 7987, 8385,
                    "1", "DENV3_29_LEFT-DENV3_29_RIGHT", 8608, 9003,
                    "1", "DENV3_31_LEFT-DENV3_31_RIGHT", 9224, 9620,
                    "1", "DENV3_33_LEFT-DENV3_33_RIGHT", 9826, 10223,
                    "1", "DENV3_35_LEFT-DENV3_35_RIGHT", 10301, 10690)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV3_2_LEFT-DENV3_2_RIGHT", 311, 711,
                    "2", "DENV3_4_LEFT-DENV3_4_RIGHT", 936, 1339,
                    "2", "DENV3_6_LEFT-DENV3_6_RIGHT", 1543, 1963,
                    "2", "DENV3_8_LEFT-DENV3_8_RIGHT", 2116, 2509,
                    "2", "DENV3_10_LEFT-DENV3_10_RIGHT", 2728, 3136,
                    "2", "DENV3_12_LEFT-DENV3_12_RIGHT", 3361, 3747,
                    "2", "DENV3_14_LEFT-DENV3_14_RIGHT", 3948, 4338,
                    "2", "DENV3_16_LEFT-DENV3_16_RIGHT", 4549, 4959,
                    "2", "DENV3_18_LEFT-DENV3_18_RIGHT", 5171, 5560,
                    "2", "DENV3_20_LEFT-DENV3_20_RIGHT", 5786, 6184,
                    "2", "DENV3_22_LEFT-DENV3_22_RIGHT", 6419, 6829,
                    "2", "DENV3_24_LEFT-DENV3_24_RIGHT", 7046, 7447,
                    "2", "DENV3_26_LEFT-DENV3_26_RIGHT", 7668, 8074,
                    "2", "DENV3_28_LEFT-DENV3_28_RIGHT", 8300, 8697,
                    "2", "DENV3_30_LEFT-DENV3_30_RIGHT", 8917, 9309,
                    "2", "DENV3_32_LEFT-DENV3_32_RIGHT", 9529, 9926,
                    "2", "DENV3_34_LEFT-DENV3_34_RIGHT", 10118, 10516)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001475.2
                         "5'UTR", 1, 94,
                         "prM", 437, 934,
                         "E", 935, 2413,
                         "NS1", 2414, 3469,
                         "NS2A", 3470, 4123,
                         "NS3", 4514, 6370,
                         "NS4A", 6371, 6751,
                         "NS4B", 6821, 7564,
                         "NS5", 7565, 10264,
                         "3'UTR", 10268, 10707)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001475.2
                         "C", 95, 436,
                         "NS2B", 4124, 4513,
                         "2K", 6752, 6820)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv3-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENGUESEQ4/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10649),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV4_1_LEFT-DENV4_1_RIGHT", 17, 508,
                    "1", "DENV4_3_LEFT-DENV4_3_RIGHT", 731, 1147,
                    "1", "DENV4_5_LEFT-DENV4_5_RIGHT", 1365, 1762,
                    "1", "DENV4_7_LEFT-DENV4_7_RIGHT", 1882, 2300,
                    "1", "DENV4_9_LEFT-DENV4_9_RIGHT", 2521, 2927,
                    "1", "DENV4_11_LEFT-DENV4_11_RIGHT", 3159, 3557,
                    "1", "DENV4_13_LEFT-DENV4_13_RIGHT", 3760, 4164,
                    "1", "DENV4_15_LEFT-DENV4_15_RIGHT", 4394, 4789,
                    "1", "DENV4_17_LEFT-DENV4_17_RIGHT", 5025, 5419,
                    "1", "DENV4_19_LEFT-DENV4_19_RIGHT", 5639, 6051,
                    "1", "DENV4_21_LEFT-DENV4_21_RIGHT", 6255, 6667,
                    "1", "DENV4_23_LEFT-DENV4_23_RIGHT", 6772, 7189,
                    "1", "DENV4_25_LEFT-DENV4_25_RIGHT", 7386, 7781,
                    "1", "DENV4_27_LEFT-DENV4_27_RIGHT", 8002, 8404,
                    "1", "DENV4_29_LEFT-DENV4_29_RIGHT", 8603, 9021,
                    "1", "DENV4_31_LEFT-DENV4_31_RIGHT", 9213, 9628,
                    "1", "DENV4_33_LEFT-DENV4_33_RIGHT", 9727, 10123,
                    "1", "DENV4_35_LEFT-DENV4_35_RIGHT", 10173, 10570)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV4_2_LEFT-DENV4_2_RIGHT", 428, 831,
                    "2", "DENV4_4_LEFT-DENV4_4_RIGHT", 1037, 1449,
                    "2", "DENV4_6_LEFT-DENV4_6_RIGHT", 1622, 2022,
                    "2", "DENV4_8_LEFT-DENV4_8_RIGHT", 2202, 2618,
                    "2", "DENV4_10_LEFT-DENV4_10_RIGHT", 2846, 3241,
                    "2", "DENV4_12_LEFT-DENV4_12_RIGHT", 3447, 3844,
                    "2", "DENV4_14_LEFT-DENV4_14_RIGHT", 4082, 4477,
                    "2", "DENV4_16_LEFT-DENV4_16_RIGHT", 4708, 5112,
                    "2", "DENV4_18_LEFT-DENV4_18_RIGHT", 5331, 5739,
                    "2", "DENV4_20_LEFT-DENV4_20_RIGHT", 5943, 6351,
                    "2", "DENV4_22_LEFT-DENV4_22_RIGHT", 6500, 6910,
                    "2", "DENV4_24_LEFT-DENV4_24_RIGHT", 7104, 7503,
                    "2", "DENV4_26_LEFT-DENV4_26_RIGHT", 7683, 8097,
                    "2", "DENV4_28_LEFT-DENV4_28_RIGHT", 8308, 8708,
                    "2", "DENV4_30_LEFT-DENV4_30_RIGHT", 8925, 9332,
                    "2", "DENV4_32_LEFT-DENV4_32_RIGHT", 9413, 9814,
                    "2", "DENV4_34_LEFT-DENV4_34_RIGHT", 10042, 10435)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1
                         "5'UTR", 1, 101,
                         "prM", 441, 938,
                         "E", 939, 2423,
                         "NS1", 2424, 3479,
                         "NS2A", 3480, 4133,
                         "NS3", 4524, 6377,
                         "NS4A", 6378, 6758,
                         "NS4B", 6828, 7562,
                         "NS5", 7563, 10262,
                         "3'UTR", 10266, 10649)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1
                         "C", 102, 440,
                         "NS2B", 4134, 4523,
                         "2K", 6759, 6827)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv4-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "ChikAsianECSA/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 11812),
                       expand = expansion(0, 0), limits = c(0, 12000)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "CHIK_400_1_LEFT_3-CHIK_400_1_RIGHT_3", 115, 451,
                    "1", "CHIK_400_3_LEFT_0-CHIK_400_3_RIGHT_0", 595, 924,
                    "1", "CHIK_400_5_LEFT_0-CHIK_400_5_RIGHT_0", 1153, 1430,
                    "1", "CHIK_400_7_LEFT_0-CHIK_400_7_RIGHT_0", 1614, 1930,
                    "1", "CHIK_400_9_LEFT_0-CHIK_400_9_RIGHT_0", 2105, 2485,
                    "1", "CHIK_400_11_LEFT_0-CHIK_400_11_RIGHT_0", 2663, 2938,
                    "1", "CHIK_400_13_LEFT_0-CHIK_400_13_RIGHT_0", 3129, 3453,
                    "1", "CHIK_400_15_LEFT_4-CHIK_400_15_RIGHT_4", 3609, 3947,
                    "1", "CHIK_400_17_LEFT_0-CHIK_400_17_RIGHT_0", 4166, 4474,
                    "1", "CHIK_400_19_LEFT_0-CHIK_400_19_RIGHT_0", 4680, 4994,
                    "1", "CHIK_400_21_LEFT_0-CHIK_400_21_RIGHT_0", 4301, 5513,
                    "1", "CHIK_400_23_LEFT_0-CHIK_400_23_RIGHT_0", 5628, 6013,
                    "1", "CHIK_400_25_LEFT_0-CHIK_400_25_RIGHT_0", 6146, 6497,
                    "1", "CHIK_400_27_LEFT_0-CHIK_400_27_RIGHT_0", 6657, 7033,
                    "1", "CHIK_400_29_LEFT_0-CHIK_400_29_RIGHT_0", 7168, 7494,
                    "1", "CHIK_400_31_LEFT_0-CHIK_400_31_RIGHT_0", 7684, 8002,
                    "1", "CHIK_400_33_LEFT_2-CHIK_400_33_RIGHT_2", 8210, 8532,
                    "1", "CHIK_400_35_LEFT_0-CHIK_400_35_RIGHT_0", 8685, 9010,
                    "1", "CHIK_400_37_LEFT_0-CHIK_400_37_RIGHT_0", 9198, 9543,
                    "1", "CHIK_400_39_LEFT_0-CHIK_400_39_RIGHT_0", 9670, 10018,
                    "1", "CHIK_400_41_LEFT_0-CHIK_400_41_RIGHT_0", 10181, 10543,
                    "1", "CHIK_400_43_LEFT_0-CHIK_400_43_RIGHT_0", 10719, 11017)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 12000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "CHIK_400_2_LEFT_0-CHIK_400_2_RIGHT_0", 294, 685,
                    "2", "CHIK_400_4_LEFT_4-CHIK_400_4_RIGHT_4", 860, 1217,
                    "2", "CHIK_400_6_LEFT_4-CHIK_400_6_RIGHT_4", 1349, 1710,
                    "2", "CHIK_400_8_LEFT_1-CHIK_400_8_RIGHT_1", 1816, 2214,
                    "2", "CHIK_400_10_LEFT_0-CHIK_400_10_RIGHT_0", 2368, 2731,
                    "2", "CHIK_400_12_LEFT_1-CHIK_400_12_RIGHT_1", 2845, 3242,
                    "2", "CHIK_400_14_LEFT_0-CHIK_400_14_RIGHT_0", 3367, 3719,
                    "2", "CHIK_400_16_LEFT_0-CHIK_400_16_RIGHT_0", 3852, 4243,
                    "2", "CHIK_400_18_LEFT_0-CHIK_400_18_RIGHT_0", 4393, 4753,
                    "2", "CHIK_400_20_LEFT_0-CHIK_400_20_RIGHT_0", 4898, 5209,
                    "2", "CHIK_400_22_LEFT_0-CHIK_400_22_RIGHT_0", 5406, 5722,
                    "2", "CHIK_400_24_LEFT_0-CHIK_400_24_RIGHT_0", 5928, 6248,
                    "2", "CHIK_400_26_LEFT_3-CHIK_400_26_RIGHT_3", 6380, 6744,
                    "2", "CHIK_400_28_LEFT_2-CHIK_400_28_RIGHT_2", 6916, 7277,
                    "2", "CHIK_400_30_LEFT_0-CHIK_400_30_RIGHT_0", 7413, 7769,
                    "2", "CHIK_400_32_LEFT_0-CHIK_400_32_RIGHT_0", 7929, 8279,
                    "2", "CHIK_400_34_LEFT_0-CHIK_400_34_RIGHT_0", 8441, 8771,
                    "2", "CHIK_400_36_LEFT_1-CHIK_400_36_RIGHT_1", 8912, 9293,
                    "2", "CHIK_400_38_LEFT_0-CHIK_400_38_RIGHT_0", 9467, 9754,
                    "2", "CHIK_400_40_LEFT_0-CHIK_400_40_RIGHT_0", 9952, 10294,
                    "2", "CHIK_400_42_LEFT_0-CHIK_400_42_RIGHT_0", 10446, 10805,
                    "2", "CHIK_400_44_LEFT_0-CHIK_400_44_RIGHT_0", 10911, 11308)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 12000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/KP164568.1
                         "5'UTR", 1, 77,
                         "NS1", 78, 1682,
                         "NS2", 1683, 4076,
                         "NS3", 4077, 5648,
                         "NS4", 5667, 7499,
                         "3'UTR", 11315, 11812)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 12000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/KP164568.1
                         "C", 7568, 8359,
                         "E3", 8360, 8542,
                         "E2", 8543, 9811,
                         "6K", 9812, 9994,
                         "E1", 9995, 11311)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 12000)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".chikv-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "CCEMHTLV1/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 9068),
                       expand = expansion(0, 0), limits = c(0, 9100)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "HTLV-1_1_LEFT-HTLV-1_1_RIGHT", 102, 505,
                    "1", "HTLV-1_3_LEFT-HTLV-1_3_RIGHT", 664, 1056,
                    "1", "HTLV-1_5_LEFT-HTLV-1_5_RIGHT", 1174, 1563,
                    "1", "HTLV-1_7_LEFT-HTLV-1_7_RIGHT", 1794, 2189,
                    "1", "HTLV-1_9_LEFT-HTLV-1_9_RIGHT", 2407, 2806,
                    "1", "HTLV-1_11_LEFT-HTLV-1_11_RIGHT", 3047, 3454,
                    "1", "HTLV-1_13_LEFT-HTLV-1_13_RIGHT", 3676, 4073,
                    "1", "HTLV-1_15_LEFT-HTLV-1_15_RIGHT", 4296, 4692,
                    "1", "HTLV-1_17_LEFT-HTLV-1_17_RIGHT", 4904, 5298,
                    "1", "HTLV-1_19_LEFT-HTLV-1_19_RIGHT", 5539, 5926,
                    "1", "HTLV-1_21_LEFT-HTLV-1_21_RIGHT", 6102, 6509,
                    "1", "HTLV-1_23_LEFT-HTLV-1_23_RIGHT", 6678, 7060,
                    "1", "HTLV-1_25_LEFT-HTLV-1_25_RIGHT", 7286, 7687,
                    "1", "HTLV-1_27_LEFT-HTLV-1_27_RIGHT", 7898, 8312,
                    "1", "HTLV-1_29_LEFT-HTLV-1_29_LEFT", 8559, 8948)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "HTLV-1_2_LEFT-HTLV-1_2_RIGHT", 349, 750,
                    "2", "HTLV-1_4_LEFT-HTLV-1_4_RIGHT", 944, 1343,
                    "2", "HTLV-1_6_LEFT-HTLV-1_6_RIGHT", 1479, 1870,
                    "2", "HTLV-1_8_LEFT-HTLV-1_8_RIGHT", 2091, 2481,
                    "2", "HTLV-1_10_LEFT-HTLV-1_10_RIGHT", 2729, 3132,
                    "2", "HTLV-1_12_LEFT-HTLV-1_12_RIGHT", 3359, 3766,
                    "2", "HTLV-1_14_LEFT-HTLV-1_14_RIGHT", 3991, 4380,
                    "2", "HTLV-1_16_LEFT-HTLV-1_16_RIGHT", 4611, 5000,
                    "2", "HTLV-1_18_LEFT-HTLV-1_18_RIGHT", 5204, 5615,
                    "2", "HTLV-1_20_LEFT-HTLV-1_20_RIGHT", 5833, 6245,
                    "2", "HTLV-1_22_LEFT-HTLV-1_22_RIGHT", 6361, 6764,
                    "2", "HTLV-1_24_LEFT-HTLV-1_24_RIGHT", 6970, 7379,
                    "2", "HTLV-1_26_LEFT-HTLV-1_26_RIGHT", 7589, 7986,
                    "2", "HTLV-1_28_LEFT-HTLV-1_28_RIGHT", 8234, 8637)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/J02029.1
                         "gag-pro-pol", 824, 5210,
                         "pX", 6857, 8382)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/J02029.1
                         "5'LTR", 23, 777,
                         "env", 5203, 6669,
                         "3'LTR", 8301, 9055)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".htlv1-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "WNV400/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11029),
                       expand = expansion(0, 0), limits = c(0, 11060)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "WNV_400_1_LEFT-WNV_400_1_RIGHT", 8, 381,
                    "1", "WNV_400_3_LEFT-WNV_400_3_RIGHT", 569, 992,
                    "1", "WNV_400_5_LEFT-WNV_400_5_RIGHT", 1209, 1616,
                    "1", "WNV_400_7_LEFT-WNV_400_7_RIGHT", 1787, 2174,
                    "1", "WNV_400_9_LEFT-WNV_400_9_RIGHT", 2332, 2739,
                    "1", "WNV_400_11_LEFT-WNV_400_11_RIGHT", 2932, 3298,
                    "1", "WNV_400_13_LEFT-WNV_400_13_RIGHT", 3449, 3836,
                    "1", "WNV_400_15_LEFT-WNV_400_15_RIGHT", 3977, 4350,
                    "1", "WNV_400_17_LEFT-WNV_400_17_RIGHT", 4564, 4962,
                    "1", "WNV_400_19_LEFT-WNV_400_19_RIGHT", 5141, 5509,
                    "1", "WNV_400_21_LEFT-WNV_400_21_RIGHT", 5697, 6093,
                    "1", "WNV_400_23_LEFT-WNV_400_23_RIGHT", 6294, 6668,
                    "1", "WNV_400_25_LEFT-WNV_400_25_RIGHT", 6881, 7295,
                    "1", "WNV_400_27_LEFT-WNV_400_27_RIGHT", 7450, 7866,
                    "1", "WNV_400_29_LEFT-WNV_400_29_RIGHT", 8040, 8449,
                    "1", "WNV_400_31_LEFT-WNV_400_31_RIGHT", 8636, 9039,
                    "1", "WNV_400_33_LEFT-WNV_400_33_RIGHT", 9233, 9595,
                    "1", "WNV_400_35_LEFT-WNV_400_35_RIGHT", 9807, 10199,
                    "1", "WNV_400_37_LEFT-WNV_400_37_RIGHT", 10375, 10803)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 11060)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "WNV_400_2_LEFT-WNV_400_2_RIGHT", 249, 680,
                    "2", "WNV_400_4_LEFT-WNV_400_4_RIGHT", 883, 1311,
                    "2", "WNV_400_6_LEFT-WNV_400_6_RIGHT", 1511, 1892,
                    "2", "WNV_400_8_LEFT-WNV_400_8_RIGHT", 2064, 2456,
                    "2", "WNV_400_10_LEFT-WNV_400_10_RIGHT", 2636, 3042,
                    "2", "WNV_400_12_LEFT-WNV_400_12_RIGHT", 3165, 3562,
                    "2", "WNV_400_14_LEFT-WNV_400_14_RIGHT", 3727, 4089,
                    "2", "WNV_400_16_LEFT-WNV_400_16_RIGHT", 4247, 4668,
                    "2", "WNV_400_18_LEFT-WNV_400_18_RIGHT", 4857, 5247,
                    "2", "WNV_400_20_LEFT-WNV_400_20_RIGHT", 5413, 5805,
                    "2", "WNV_400_22_LEFT-WNV_400_22_RIGHT", 5996, 6394,
                    "2", "WNV_400_24_LEFT-WNV_400_24_RIGHT", 6563, 6988,
                    "2", "WNV_400_26_LEFT-WNV_400_26_RIGHT", 7175, 7555,
                    "2", "WNV_400_28_LEFT-WNV_400_28_RIGHT", 7757, 8137,
                    "2", "WNV_400_30_LEFT-WNV_400_30_RIGHT", 8349, 8731,
                    "2", "WNV_400_32_LEFT-WNV_400_32_RIGHT", 8928, 9329,
                    "2", "WNV_400_34_LEFT-WNV_400_34_RIGHT", 9500, 9920,
                    "2", "WNV_400_36_LEFT-WNV_400_36_RIGHT", 10086, 10514,
                    "2", "WNV_400_38_LEFT-WNV_400_38_RIGHT", 10591, 10969)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 11060)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1
                         "5'UTR", 1, 96,
                         "pr", 466, 741,
                         "M", 742, 966,
                         "E", 967, 2469,
                         "NS1", 2470, 3525,
                         "NS2A", 3526, 4218,
                         "NS3", 4612, 6468,
                         "NS4A", 6469, 6846,
                         "NS4B", 6916, 7680,
                         "NS5", 7681, 10395,
                         "3'UTR", 10399, 11029)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 11060)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1
                         "C", 97, 411,
                         "NS2B", 4219, 4611,
                         "2K", 6847, 6915)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 11060)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".wnv-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "HIV1Sanabani2006/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 2000, 4000, 6000, 8000, 9719),
                       expand = expansion(0, 0), limits = c(0, 9800)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "SC-ANS_F-SC-ANA_R", 545, 2610,
                    "1", "IN-B01S_F-C-BNA_R", 3235, 5220,
                    "1", "SC-DNS_F-SC-DNA_R", 7718, 9537)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "SC-BNS_F-IN-B01R_R", 2156, 3791,
                    "2", "SC-CNS_F-SC-CNA_R", 4889, 7808)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html
                         "5'LTR", 1, 634,
                         "gag", 790, 2292,
                         "vif", 5041, 5619,
                         "tat", 8379, 8424,
                         "nef", 8797, 9417)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html
                         "tat", 5831, 6045,
                         "vpu", 6062, 6310,
                         "rev", 8379, 8653,
                         "3'LTR", 9086, 9719)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome3 <- tribble(~"gene", ~"start", ~"end", # https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html
                         "pol", 2085, 5096,
                         "vpr", 5559, 5850,
                         "rev", 5970, 6045,
                         "env", 6225, 8795)
  map2plot3 <- map2genome3 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9800)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")  
  output <-  paste0(output, ".hiv1-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 / map2plot3 + plot_layout(nrow = 7, heights = c(3, .1, .1, .1, .3, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "RSVA/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000),
                       expand = expansion(0, 0), limits = c(0, 15400)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "A1f-A1r", 1, 1834,
                    "1", "AB3f-A3r", 2926, 4885,
                    "1", "A5f-A5r", 6394, 7976,
                    "1", "A7f-A7r", 9450, 11045,
                    "1", "A9f-A9r", 12343, 14144)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "A2f-A2r", 1584, 3459,
                    "2", "A4f-A4r", 4732, 6586,
                    "2", "A6f-A6r", 7670, 9569,
                    "2", "A8f-A8r", 10731, 12668,
                    "2", "A10f-A10r", 13819, 15277)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/MZ515716.1
                        "NS1", 99, 518,
                        "NS2", 626, 1000,
                        "N", 1139, 2314,
                        "P", 2347, 3072,
                        "M", 3262, 4032,
                        "SH", 4301, 4498,
                        "G", 4689, 5642,
                        "F", 5719, 7443,
                        "M2", 7670, 8495,
                        "L", 8561, 15061)
  map2plot <- map2genome %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".rsva-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "RSVA/V2") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000),
                       expand = expansion(0, 0), limits = c(0, 15400)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "A1f-A1r", 1, 1834,
                    "1", "AB3f-A3r", 2926, 4885,
                    "1", "A5f-A5r", 6394, 7976,
                    "1", "A7f-A7r", 9450, 11045,
                    "1", "A9f-A9r", 12343, 14144)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "A2f-A2r", 1584, 3459,
                    "2", "A4f-A4r", 4732, 6586,
                    "2", "A6f-A6r", 7670, 9569,
                    "2", "A8f-A8r", 10731, 12668,
                    "2", "A10f-A10r", 13819, 15277)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/MZ515716.1
                        "NS1", 99, 518,
                        "NS2", 626, 1000,
                        "N", 1139, 2314,
                        "P", 2347, 3072,
                        "M", 3262, 4032,
                        "SH", 4301, 4498,
                        "G", 4689, 5642,
                        "F", 5719, 7443,
                        "M2", 7670, 8495,
                        "L", 8561, 15061)
  map2plot <- map2genome %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".rsva-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "RSVA/V3") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000),
                       expand = expansion(0, 0), limits = c(0, 15400)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "A1f-A1r", 1, 1834,
                    "1", "AB3f-A3r", 2926, 4885,
                    "1", "A5f-A5r", 6394, 7976,
                    "1", "A7f-A7r", 9450, 11045,
                    "1", "A9f-A9r", 12343, 14144)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "A2f-A2r", 1584, 3459,
                    "2", "A4f-A4r", 4732, 6586,
                    "2", "A6f-A6r", 7670, 9569,
                    "2", "A8f-A8r", 10731, 12668,
                    "2", "A10f2-A10r", 13680, 15277)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/MZ515716.1
                        "NS1", 99, 518,
                        "NS2", 626, 1000,
                        "N", 1139, 2314,
                        "P", 2347, 3072,
                        "M", 3262, 4032,
                        "SH", 4301, 4498,
                        "G", 4689, 5642,
                        "F", 5719, 7443,
                        "M2", 7670, 8495,
                        "L", 8561, 15061)
  map2plot <- map2genome %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".rsva-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "RSVB/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000),
                       expand = expansion(0, 0), limits = c(0, 15400)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "B1f-B1r", 1, 1704,
                    "1", "AB3f-B3r", 2926, 4352,
                    "1", "B5f-B5r", 5743, 7421,
                    "1", "B7f-B7r", 8789, 10551,
                    "1", "B9f-B9r", 11903, 13694)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "B2f-B2r", 1354, 3025,
                    "2", "B4f-B4r", 3985, 5832,
                    "2", "B6f-B6r", 7156, 8878,
                    "2", "B8f-B8r", 10274, 12234,
                    "2", "B10f-B10r", 13332, 15278)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/MN163126.1
                        "NS1", 100, 519,
                        "NS2", 627, 1001,
                        "N", 1139, 2314,
                        "P", 2347, 3072,
                        "M", 3262, 4032,
                        "SH", 4301, 4498,
                        "G", 4687, 5619,
                        "F", 5717, 7441,
                        "M2", 7668, 8493,
                        "L", 8559, 15059)
  map2plot <- map2genome %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".rsvb-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "RSVB/V2") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000),
                       expand = expansion(0, 0), limits = c(0, 15400)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "B1f-B1r", 1, 1704,
                    "1", "AB3f-B3r", 2926, 4352,
                    "1", "B5f-B5r", 5743, 7421,
                    "1", "B7f-B7r", 8789, 10551,
                    "1", "B9f-B9r", 11903, 13694)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "B2f-B2r", 1354, 3025,
                    "2", "B4f-B4r", 3985, 5832,
                    "2", "B6f-B6r", 7156, 8878,
                    "2", "B8f-B8r", 10274, 12234,
                    "2", "B10f-B10r", 13332, 15278)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/MN163126.1
                        "NS1", 100, 519,
                        "NS2", 627, 1001,
                        "N", 1139, 2314,
                        "P", 2347, 3072,
                        "M", 3262, 4032,
                        "SH", 4301, 4498,
                        "G", 4687, 5619,
                        "F", 5717, 7441,
                        "M2", 7668, 8493,
                        "L", 8559, 15059)
  map2plot <- map2genome %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".rsvb-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "RSVB/V3") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000),
                       expand = expansion(0, 0), limits = c(0, 15400)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "B1f-B1r", 1, 1704,
                    "1", "AB3f-B3r", 2926, 4352,
                    "1", "B5f-B5r", 5743, 7421,
                    "1", "B7f-B7r", 8789, 10551,
                    "1", "B9f-B9r", 11903, 13694)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "B2f-B2r", 1354, 3025,
                    "2", "B4f-B4r", 3985, 5832,
                    "2", "B6f-B6r", 7156, 8878,
                    "2", "B8f-B8r", 10274, 12234,
                    "2", "B10f-B10r", 13332, 15278)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/MN163126.1
                        "NS1", 100, 519,
                        "NS2", 627, 1001,
                        "N", 1139, 2314,
                        "P", 2347, 3072,
                        "M", 3262, 4032,
                        "SH", 4301, 4498,
                        "G", 4687, 5619,
                        "F", 5717, 7441,
                        "M2", 7668, 8493,
                        "L", 8559, 15059)
  map2plot <- map2genome %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 15400)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".rsvb-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "HTLV1DemincoF/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 9068),
                       expand = expansion(0, 0), limits = c(0, 9100)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000),
                       expand = expansion(0, 0), limits = c(1, 2001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
    map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "AV1 Fn1a+/b+ - AV1 R", 23, 4748)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "BV1 F2 - BV1 Ra+/bn+", 4124, 9056)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 3) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/J02029.1
                         "gag-pro-pol", 824, 5210,
                         "pX", 6857, 8382)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/J02029.1
                         "5'LTR", 23, 777,
                         "env", 5203, 6669,
                         "3'LTR", 8301, 9055)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 9100)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".htlv1-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENV1CADDE/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10735),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV1_1_LEFT-DENV1_1_RIGHT", 134, 1080,
                    "1", "DENV1_3_LEFT-DENV1_3_RIGHT", 1447, 2380,
                    "1", "DENV1_5_LEFT-DENV1_5_RIGHT", 2867, 3778,
                    "1", "DENV1_7_LEFT-DENV1_7_RIGHT", 4130, 5011,
                    "1", "DENV1_9_LEFT-DENV1_9_RIGHT", 5402, 6252,
                    "1", "DENV1_11_LEFT-DENV1_11_RIGHT", 6784, 7701,
                    "1", "DENV1_13_LEFT-DENV1_13_RIGHT", 8151, 8990,
                    "1", "DENV1_15_LEFT-DENV1_15_RIGHT", 9353, 10182)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV1_2_LEFT-DENV1_2_RIGHT", 825, 1648,
                    "2", "DENV1_4_LEFT-DENV1_4_RIGHT", 2110, 2986,
                    "2", "DENV1_6_LEFT-DENV1_6_RIGHT", 3279, 4242,
                    "2", "DENV1_8_LEFT-DENV1_8_RIGHT", 4749, 5565,
                    "2", "DENV1_10_LEFT-DENV1_10_RIGHT", 6064, 6953,
                    "2", "DENV1_12_LEFT-DENV1_12_RIGHT", 7590, 8519,
                    "2", "DENV1_14_LEFT-DENV1_14_RIGHT", 8829, 9755)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1
                         "5'UTR", 1, 94,
                         "prM", 437, 934,
                         "E", 935, 2419,
                         "NS1", 2420, 3475,
                         "NS2A", 3476, 4129,
                         "NS3", 4520, 6376,
                         "NS4A", 6377, 6757,
                         "NS4B", 6827, 7573,
                         "NS5", 7574, 10270,
                         "3'UTR", 10274, 10735)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1
                         "C", 95, 436,
                         "NS2B", 4130, 4519,
                         "2K", 6758, 6826)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv1-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENV2CADDE/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV2_1_LEFT2-DENV2_1_RIGHT", 85, 1046,
                    "1", "DENV2_3_LEFT-DENV2_3_RIGHT", 1642, 2542,
                    "1", "DENV2_5_LEFT-DENV2_5_RIGHT", 2708, 3760,
                    "1", "DENV2_7_LEFT-DENV2_7_RIGHT", 4345, 5264,
                    "1", "DENV2_9_LEFT-DENV2_9_RIGHT", 5860, 6812,
                    "1", "DENV2_11_LEFT-DENV2_11_RIGHT", 7519, 8572,
                    "1", "DENV2_13_LEFT-DENV2_13_RIGHT", 9181, 10113)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV2_2_LEFT-DENV2_2_RIGHT", 920, 1805,
                    "2", "DENV2_4_LEFT-DENV2_4_RIGHT", 2423, 3363,
                    "2", "DENV2_6_LEFT-DENV2_6_RIGHT", 3419, 4421,
                    "2", "DENV2_8_LEFT-DENV2_8_RIGHT", 5051, 5949,
                    "2", "DENV2_10_LEFT-DENV2_10_RIGHT", 6741, 7679,
                    "2", "DENV2_12_LEFT-DENV2_12_RIGHT", 8427, 9333)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001474.2
                         "5'UTR", 1, 96,
                         "prM", 439, 936,
                         "E", 937, 2421,
                         "NS1", 2422, 3477,
                         "NS2A", 3478, 4131,
                         "NS3", 4522, 6375,
                         "NS4A", 6376, 6756,
                         "NS4B", 6826, 7569,
                         "NS5", 7570, 10269,
                         "3'UTR", 10273, 10723)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001474.2
                         "C", 97, 438,
                         "NS2B", 4132, 4521,
                         "2K", 6757, 6825)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv2-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENV3CADDE/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV3_1_LEFT-DENV3_1_RIGHT", 132, 1084,
                    "1", "DENV3_3_LEFT-DENV3_3_RIGHT", 1303, 2142,
                    "1", "DENV3_5_LEFT-DENV3_5_RIGHT", 2773, 3623,
                    "1", "DENV3_7_LEFT-DENV3_7_RIGHT", 4011, 4825,
                    "1", "DENV3_9_LEFT-DENV3_9_RIGHT", 5252, 6135,
                    "1", "DENV3_11_LEFT-DENV3_11_RIGHT", 6605, 7449,
                    "1", "DENV3_13_LEFT-DENV3_13_RIGHT", 8029, 8887,
                    "1", "DENV3_15_LEFT-DENV3_15_RIGHT", 9330, 10140)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV3_2_LEFT-DENV3_2_RIGHT", 946, 1892,
                    "2", "DENV3_4_LEFT-DENV3_4_RIGHT", 2033, 2882,
                    "2", "DENV3_6_LEFT-DENV3_6_RIGHT", 3480, 4322,
                    "2", "DENV3_8_LEFT-DENV3_8_RIGHT", 4701, 5559,
                    "2", "DENV3_10_LEFT-DENV3_10_RIGHT", 5965, 6778,
                    "2", "DENV3_12_LEFT-DENV3_12_RIGHT", 7259, 8142,
                    "2", "DENV3_14_LEFT-DENV3_14_RIGHT", 8706, 9654)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001475.2
                         "5'UTR", 1, 94,
                         "prM", 437, 934,
                         "E", 935, 2413,
                         "NS1", 2414, 3469,
                         "NS2A", 3470, 4123,
                         "NS3", 4514, 6370,
                         "NS4A", 6371, 6751,
                         "NS4B", 6821, 7564,
                         "NS5", 7565, 10264,
                         "3'UTR", 10268, 10707)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_001475.2
                         "C", 95, 436,
                         "NS2B", 4124, 4513,
                         "2K", 6752, 6820)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv3-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENV4CADDE/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10649),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV4_1_LEFT-DENV4_1_RIGHT", 40, 787,
                    "1", "DENV4_3_LEFT-DENV4_3_RIGHT", 1272, 2082,
                    "1", "DENV4_5_LEFT-DENV4_5_RIGHT", 2496, 3332,
                    "1", "DENV4_7_LEFT-DENV4_7_RIGHT", 3536, 4398,
                    "1", "DENV4_9_LEFT-DENV4_9_RIGHT", 4753, 5618,
                    "1", "DENV4_11_LEFT-DENV4_11_RIGHT", 6042, 6766,
                    "1", "DENV4_13_LEFT-DENV4_13_RIGHT", 7277, 8091,
                    "1", "DENV4_15_LEFT-DENV4_15_RIGHT", 8487, 9312)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV4_2_LEFT-DENV4_2_RIGHT", 622, 1386,
                    "2", "DENV4_4_LEFT-DENV4_4_RIGHT", 1870, 2620,
                    "2", "DENV4_6_LEFT-DENV4_6_RIGHT", 3081, 3812,
                    "2", "DENV4_8_LEFT-DENV4_8_RIGHT", 4226, 5066,
                    "2", "DENV4_10_LEFT-DENV4_10_RIGHT", 5509, 6260,
                    "2", "DENV4_12_LEFT-DENV4_12_RIGHT", 6605, 7422,
                    "2", "DENV4_14_LEFT-DENV4_14_RIGHT", 7929, 8758,
                    "2", "DENV4_16_LEFT-DENV4_16_RIGHT", 9071, 9912)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1
                         "5'UTR", 1, 101,
                         "prM", 441, 938,
                         "E", 939, 2423,
                         "NS1", 2424, 3479,
                         "NS2A", 3480, 4133,
                         "NS3", 4524, 6377,
                         "NS4A", 6378, 6758,
                         "NS4B", 6828, 7562,
                         "NS5", 7563, 10262,
                         "3'UTR", 10266, 10649)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1
                         "C", 102, 440,
                         "NS2B", 4134, 4523,
                         "2K", 6759, 6827)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv4-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "OROVFN400L/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 2000, 3000, 4000, 5000, 6846),
                       expand = expansion(0, 0), limits = c(0, 6870)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "OROV_L_v2_0_LEFT-OROV_L_v2_0_RIGHT", 9, 400,
                    "2", "OROV_L_v2_2_LEFT-OROV_L_v2_2_RIGHT", 335, 738,
                    "2", "OROV_L_v2_4_LEFT-OROV_L_v2_4_RIGHT", 943, 1352,
                    "2", "OROV_L_v2_6_LEFT-OROV_L_v2_6_RIGHT", 1536, 1933,
                    "2", "OROV_L_v2_8_LEFT-OROV_L_v2_8_RIGHT", 2129, 2542,
                    "2", "OROV_L_v2_10_LEFT/alt/alt1-OROV_L_v2_10_RIGHT/alt1", 2751, 3161,
                    "2", "OROV_L_v2_12_LEFT/alt-OROV_L_v2_12_RIGHT", 3380, 3781,
                    "2", "OROV_L_v2_14_LEFT-OROV_L_v2_14_RIGHT", 3987, 4383,
                    "2", "OROV_L_v2_16_LEFT-OROV_L_v2_16_RIGHT", 4571, 4982,
                    "2", "OROV_L_v2_18_LEFT-OROV_L_v2_18_RIGHT/alt", 5172, 5585,
                    "2", "OROV_L_v2_20_LEFT-OROV_L_v2_20_RIGHT", 5814, 6221,
                    "2", "OROV_L_v2_22_LEFT-OROV_L_v2_22_RIGHT", 6389, 6789)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 6870)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "OROV_L_v2_1_LEFT-OROV_L_v2_1_RIGHT", 33, 440,
                    "1", "OROV_L_v2_3_LEFT/alt-OROV_L_v2_3_RIGHT", 629, 1042,
                    "1", "OROV_L_v2_5_LEFT-OROV_L_v2_5_RIGHT", 1238, 1642,
                    "1", "OROV_L_v2_7_LEFT-OROV_L_v2_7_RIGHT", 1825, 2228,
                    "1", "OROV_L_v2_9_LEFT-OROV_L_v2_9_RIGHT", 2441, 2851,
                    "1", "OROV_L_v2_11_LEFT-OROV_L_v2_11_RIGHT", 3066, 3482,
                    "1", "OROV_L_v2_13_LEFT-OROV_L_v2_13_RIGHT", 3665, 4069,
                    "1", "OROV_L_v2_15_LEFT-OROV_L_v2_15_RIGHT", 4283, 4693,
                    "1", "OROV_L_v2_17_LEFT-OROV_L_v2_17_RIGHT", 4894, 5284,
                    "1", "OROV_L_v2_19_LEFT-OROV_L_v2_19_RIGHT", 5501, 5909,
                    "1", "OROV_L_v2_21_LEFT-OROV_L_v2_21_RIGHT", 6103, 6505,
                    "1", "OROV_L_v2_23_LEFT-OROV_L_v2_23_RIGHT", 6526, 6849)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 6870)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/KP691612.1
                         "OROVsLgp1", 44, 6802)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 6870)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".orovL-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "OROVFN400M/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4385),
                       expand = expansion(0, 0), limits = c(0, 4410)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "OROV_M_400_1_LEFT/alt1-OROV_M_400_1_RIGHT", 12, 385,
                    "1", "OROV_M_400_3_LEFT-OROV_M_400_3_RIGHT", 542, 900,
                    "1", "OROV_M_400_5_LEFT-OROV_M_400_5_RIGHT", 1111, 1490,
                    "1", "OROV_M_400_7_LEFT-OROV_M_400_7_RIGHT", 1671, 2047,
                    "1", "OROV_M_400_9_LEFT-OROV_M_400_9_RIGHT", 2178, 2578,
                    "1", "OROV_M_400_11_LEFT-OROV_M_400_11_RIGHT", 2712, 3123,
                    "1", "OROV_M_400_13_LEFT-OROV_M_400_13_RIGHT", 3302, 3699,
                    "1", "OROV_M_400_15_LEFT-OROV_M_400_15_RIGHT", 3847, 4232)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 4410)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "OROV_M_400_2_LEFT/alt1/alt2-OROV_M_400_2_RIGHT", 284, 648,
                    "2", "OROV_M_400_4_LEFT-OROV_M_400_4_RIGHT/alt1", 801, 1212,
                    "2", "OROV_M_400_6_LEFT-OROV_M_400_6_RIGHT", 1382, 1784,
                    "2", "OROV_M_400_8_LEFT/alt1-OROV_M_400_8_RIGHT/alt1/alt2", 1940, 2310,
                    "2", "OROV_M_400_10_LEFT-OROV_M_400_10_RIGHT", 2475, 2844,
                    "2", "OROV_M_400_12_LEFT-OROV_M_400_12_RIGHT", 3025, 3412,
                    "2", "OROV_M_400_14_LEFT-OROV_M_400_14_RIGHT", 3592, 3982,
                    "2", "OROV_M_400_16_LEFT-OROV_M_400_16_RIGHT", 4008, 4381)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 4410)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/KP691622.1
                         "OROVsMgp1", 32, 4294)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 4410)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".orovM-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 + plot_layout(nrow = 5, heights = c(3, .1, .1, .1, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "OROVFN400S/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 946),
                       expand = expansion(0, 0), limits = c(0, 950)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "OROV_S_400_0_LEFT-OROV_S_400_0_RIGHT", 5, 264,
                    "2", "OROV_S_400_2_LEFT-OROV_S_400_2_RIGHT", 279, 710)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 950)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "OROV_S_400_1_LEFT-OROV_S_400_1_RIGHT", 38, 398,
                    "1", "OROV_S_400_3_LEFT-OROV_S_400_3_RIGHT", 548, 909)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 950)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/KP691623.1
                         "OROVsNgp1", 45, 740)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 950)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/KP691623.1
                         "OROVsNgp2", 64,342)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 950)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
    output <-  paste0(output, ".orovS-coverage.pdf")
    plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}

if (primer_scheme == "DENV2GII2022FN/V1") {
  depcov <- ggplot() +
    geom_line(data = depth_coverage, aes(x = position, y = depth), linewidth = .4, colour = "black") +
    labs(title = paste0(id_sample), subtitle = paste0(primer_scheme_2),
         y = "Per base coverage (x)", x = NULL) +
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000, 7000, 9000, 10723),
                       expand = expansion(0, 0), limits = c(0, 10750)) +
    scale_y_continuous(breaks = c(1, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                       expand = expansion(0, 0), limits = c(1, 100001),
                       labels = function(x) format(x, scientific = FALSE), trans = "log10") +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.y = element_text(angle = 90, size = 14),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(hjust = 1, size = 9)) +
    geom_hline(yintercept = 10, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 20, linetype = "dotted", colour = "black")
  map1p1 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "1", "DENV2_GII_2022_1_LEFT-DENV2_GII_2022_1_RIGHT/alt", 4, 410,
                    "1", "DENV2_GII_2022_3_LEFT-DENV2_GII_2022_3_RIGHT", 632, 1035,
                    "1", "DENV2_GII_2022_5_LEFT-DENV2_GII_2022_5_RIGHT", 1249, 1654,
                    "1", "DENV2_GII_2022_7_LEFT-DENV2_GII_2022_7_RIGHT", 1866, 2273,
                    "1", "DENV2_GII_2022_9_LEFT-DENV2_GII_2022_9_RIGHT", 2486, 2904,
                    "1", "DENV2_GII_2022_11_LEFT-DENV2_GII_2022_11_RIGHT", 3126, 3538,
                    "1", "DENV2_GII_2022_13_LEFT-DENV2_GII_2022_13_RIGHT", 3766, 4170,
                    "1", "DENV2_GII_2022_15_LEFT-DENV2_GII_2022_15_RIGHT", 4399, 4813,
                    "1", "DENV2_GII_2022_17_LEFT-DENV2_GII_2022_17_RIGHT", 5033, 5445,
                    "1", "DENV2_GII_2022_19_LEFT-DENV2_GII_2022_19_RIGHT", 5671, 6083,
                    "1", "DENV2_GII_2022_21_LEFT-DENV2_GII_2022_21_RIGHT", 6322, 6717,
                    "1", "DENV2_GII_2022_23_LEFT-DENV2_GII_2022_23_RIGHT", 6933, 7346,
                    "1", "DENV2_GII_2022_25_LEFT-DENV2_GII_2022_25_RIGHT", 7561, 7965,
                    "1", "DENV2_GII_2022_27_LEFT-DENV2_GII_2022_27_RIGHT", 8186, 8592,
                    "1", "DENV2_GII_2022_29_LEFT-DENV2_GII_2022_29_RIGHT", 8815, 9222,
                    "1", "DENV2_GII_2022_31_LEFT-DENV2_GII_2022_31_RIGHT", 9447, 9859,
                    "1", "DENV2_GII_2022_33_LEFT-DENV2_GII_2022_33_RIGHT", 10084, 10482)
  map1plot1 <- map1p1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("1" = "red"))
  map1p2 <- tribble(~"pool", ~"amplicon", ~"start", ~"end",
                    "2", "DENV2_GII_2022_2_LEFT-DENV2_GII_2022_2_RIGHT", 314, 729,
                    "2", "DENV2_GII_2022_4_LEFT-DENV2_GII_2022_4_RIGHT", 949, 1345,
                    "2", "DENV2_GII_2022_6_LEFT-DENV2_GII_2022_6_RIGHT", 1561, 1950,
                    "2", "DENV2_GII_2022_8_LEFT-DENV2_GII_2022_8_RIGHT", 2183, 2582,
                    "2", "DENV2_GII_2022_10_LEFT-DENV2_GII_2022_10_RIGHT", 2813, 3222,
                    "2", "DENV2_GII_2022_12_LEFT-DENV2_GII_2022_12_RIGHT", 3455, 3860,
                    "2", "DENV2_GII_2022_14_LEFT-DENV2_GII_2022_14_RIGHT", 4074, 4490,
                    "2", "DENV2_GII_2022_16_LEFT-DENV2_GII_2022_16_RIGHT", 4730, 5123,
                    "2", "DENV2_GII_2022_18_LEFT-DENV2_GII_2022_18_RIGHT", 5357, 5766,
                    "2", "DENV2_GII_2022_20_LEFT-DENV2_GII_2022_20_RIGHT", 6000, 6407,
                    "2", "DENV2_GII_2022_22_LEFT-DENV2_GII_2022_22_RIGHT", 6624, 7026,
                    "2", "DENV2_GII_2022_24_LEFT-DENV2_GII_2022_24_RIGHT", 7241, 7647,
                    "2", "DENV2_GII_2022_26_LEFT-DENV2_GII_2022_26_RIGHT", 7877, 8280,
                    "2", "DENV2_GII_2022_28_LEFT-DENV2_GII_2022_28_RIGHT", 8504, 8915,
                    "2", "DENV2_GII_2022_30_LEFT-DENV2_GII_2022_30_RIGHT", 9132, 9540,
                    "2", "DENV2_GII_2022_32_LEFT-DENV2_GII_2022_32_RIGHT", 9770, 10183,
                    "2", "DENV2_GII_2022_34_LEFT-DENV2_GII_2022_34_RIGHT", 10192, 10602)
  map1plot2 <- map1p2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10, fill = pool), alpha = .4) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = amplicon), size = 1) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("2" = "blue"))
  map2genome1 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/PV789655.1
                         "5'UTR", 1, 95,
                         "prM", 439, 936,
                         "E", 937, 2421,
                         "NS1", 2422, 3477,
                         "NS2A", 3478, 4131,
                         "NS3", 4522, 6375,
                         "NS4A", 6376, 6756,
                         "NS4B", 6826, 7569,
                         "NS5", 7570, 10269,
                         "3'UTR", 10276, 10725)
  map2plot1 <- map2genome1 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  map2genome2 <- tribble(~"gene", ~"start", ~"end", # https://www.ncbi.nlm.nih.gov/nuccore/PV789655.1
                         "C", 97, 396,
                         "NS2B", 4132, 4521,
                         "2K", 6757, 6825)
  map2plot2 <- map2genome2 %>% ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 8, ymax = 10),
              linewidth = .2, fill = "green", colour = "darkgray", alpha = .3) +
    geom_text(aes(x = (start + end) / 2, y = 9, label = gene), size = 4) +
    scale_x_continuous(expand = expansion(0, 0), limits = c(0, 10750)) +
    theme_void() + theme(legend.position = "none") + coord_cartesian(clip = "off")
  output <-  paste0(output, ".denv2gii-coverage.pdf")
  plot <- depcov / map1plot1 / map1plot2 / plot_spacer() / map2plot1 / map2plot2 + plot_layout(nrow = 6, heights = c(3, .1, .1, .1, .3, .3))
  save_plot(output, plot, base_height = 7, base_width = 20)
}