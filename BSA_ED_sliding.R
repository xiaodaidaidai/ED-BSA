library(argparser, quietly = TRUE)

# Create a parser
p <- arg_parser("Mapping the distribution of ED correlation values on chromosomes")

# Add command line arguments
p <- add_argument(p, "ED", help = "Input: ED.tsv", type = "character")
p <- add_argument(p, "window", help = "Input: sliding_window.tsv", type = "character")
p <- add_argument(p, "chr", help = "Input: Chromosome list file (chr.len)", type = "character")
p <- add_argument(p, "threshold", help = "Input: Threshold of confidence interval (99 or 95)", type = "character")
p <- add_argument(p, "outpre", help = "Output prefix", type = "character")

# Parse command line arguments
argv <- parse_args(p)

# Load required libraries
library(ggplot2)
library(tidyverse)
options(scipen = 100)

fsnp <- argv$ED
fwindow <- argv$window
fchr <- argv$chr
cut <- as.numeric(argv$threshold)
outpre <- argv$outpre

snp <- read.delim(fsnp, skip = 1, header = F)
window <- read.delim(fwindow, skip = 1, header = F)
chr <- read.delim(fchr, skip = 1, header = F)

window <- window %>% as_tibble() %>% filter(V1 %in% chr$V1, !is.na(V6))
window$V1 <- factor(window$V1, levels = chr$V1)

cutvalue <- quantile(window$V6, cut/100)
print(cutvalue)

outwindow <- window %>% filter(V6 > cutvalue)
write.table(outwindow, file = paste(outpre, "greater_than_threshold_EDsliding.tsv", sep = "."), col.names = F, row.names = F, quote = F, sep = "\t")

snp <- snp %>% as_tibble() %>% filter(V1 %in% chr$V1, !is.na(V5))
snp$V1 <- factor(snp$V1, levels = chr$V1)

chr_unit <- 1000000
snp$V2 <- snp$V2/chr_unit
window$V2 <- window$V2/chr_unit

colorlist <- rep(c("Darkgreen", "DarkOrange"), times = 1, len = nrow(chr))

pbulk <- ggplot(data = snp, aes(x = V2, y = V5)) +
  geom_point(aes(color = V1), size = 0.5) +
  geom_line(data = window, aes(x = V2, y = V6), color = "black", linewidth = 0.8) +
  geom_hline(yintercept = cutvalue, color = "red", linetype = 2, linewidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_color_manual(values = colorlist) +
  facet_wrap(~V1, nrow = 1, scales = "free_x", strip.position = 'bottom') +
  ylab("Euclidean distance") +
  xlab("Chromosome") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0, "mm"),
    panel.border = element_rect(linetype = 2, fill = NA, color = "grey"),
    strip.placement = "outside",
    strip.text.x = element_text(angle = 45, hjust = 0.5),
    strip.background.x = element_blank()
  )

ggsave(paste(outpre, "plot-EDsliding.pdf", sep = "."), plot = pbulk, device = "pdf", width = 10, height = 5 )
ggsave(paste(outpre, "plot-EDsliding.png", sep = "."), plot = pbulk, device = "png", width = 10, height = 5 )
