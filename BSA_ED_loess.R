library(argparser, quietly = TRUE)

# Create a parser
p <- arg_parser("Generate SNP-index figure for mutmap.")

# Add command line arguments
p <- add_argument(p, "ED", help = "Input: ED.tsv", type = "character")
p <- add_argument(p, "chr", help = "Input: Chromosome length file (chr.len)", type = "character")
p <- add_argument(p, "threshold", help = "Input: Threshold of confidence interval (99 or 95)", type = "character")
p <- add_argument(p, "outpre", help = "Output prefix", type = "character")

# Parse command line arguments
argv <- parse_args(p)

# Load required libraries
library(ggplot2)
library(tidyverse)
library(MMAPPR2)
library(GenomicRanges)

# Read input files
fsnp <- argv$ED
fchr <- argv$chr
cut <- as.numeric(argv$threshold)
outpre <- argv$outpre

options(scipen = 100)

# Read SNP and chromosome length data
snp <- read.delim(fsnp, skip = 1, header = F)
chr <- read.delim(fchr, skip = 1, header = F)

# Filter SNP data based on chromosomes, remove NA values, and modify chromosome levels order
snp <- snp %>% 
  as_tibble() %>% 
  filter(V1 %in% chr$V1, !is.na(V5)) %>%
  mutate(V1 = factor(V1, levels = chr$V1)) %>%
  mutate(V1 = as.character(V1))

# Define classes for MMAPPR2 data
setClass("param", representation(loessOptResolution = "numeric", loessOptCutFactor = "numeric", peakIntervalWidth = "numeric"))
setClass("EDdata", representation(param = "param", distance = "list", peaks = "list", candidates = "list"))
chrRanges <- GRanges(seqnames = as.character(chr$V1), ranges = IRanges(1, chr$V2), strand = rep("+", length(chr$V1)))

chrList <- list()

# Store range for each chromosome as a list item
for (i in suppressWarnings(GenomeInfoDb::orderSeqlevels(as.character(GenomeInfoDb::seqnames(chrRanges))))) {
  chrList[[toString(GenomeInfoDb::seqnames(chrRanges[i]))]] <- chrRanges[i]
}

# Function to generate distance data
storeED <- function(chrRange, snp) {
  ED <- snp %>% dplyr::filter(V1 == as.character(GenomeInfoDb::seqnames(chrRange)))
  distanceDf <- data.frame(pos = ED$V2, distance = ED$V5, chr = ED$V1)
  resultList <- list(distanceDf = distanceDf)
  resultList$seqname <- as.character(GenomeInfoDb::seqnames(chrRange))
  return(resultList)
}

# Create MMAPPR2 data
param <- new("param", loessOptResolution = 0.01, loessOptCutFactor = 0.1, peakIntervalWidth = 0.95)
EDdata <- new("EDdata", param = param)

EDdata@distance <- BiocParallel::bplapply(chrList, storeED, snp = snp)

# Perform loess fitting with span values determined by AICc optimization
postLoessMD <- loessFit(EDdata)

# Output loess fitted table and best span values
LoessOut <- data.frame()
span <- data.frame()

for (i in names(postLoessMD@distance)) {
  CHR <- i
  dloess <- data.frame(
    CHR = rep(CHR, times = postLoessMD@distance[[i]]$loess$n),
    pos = postLoessMD@distance[[i]]$loess$x,
    ED = as.vector(postLoessMD@distance[[i]]$loess$y),
    LOESS = postLoessMD@distance[[i]]$loess$fitted
  )
  LoessOut <- rbind(LoessOut, dloess)
  s <- data.frame(CHR = CHR, SPAN = postLoessMD@distance[[i]]$loess$s)
  span <- rbind(span, s)
}

write.table(LoessOut, file = paste(outpre, "EDloessAICc.fitted.tsv", sep = "."), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(span, file = paste(outpre, "EDloessAICc.span.tsv", sep = "."), sep = "\t", quote = F, row.names = F, col.names = T)

# Read loess fitted data
filename <- paste(outpre, "EDloessAICc.fitted.tsv", sep = ".")
LoessOut <- read.delim(filename, skip = 1, header = F)

# Set chromosome unit and calculate cut value
chr_unit <- 1000000
LoessOut$V2 <- LoessOut$V2 / chr_unit
cutvalue <- quantile(LoessOut$V4, cut / 100)

# Print cut value
print(cutvalue)

# Set color list
colorlist <- rep(c("Darkgreen", "DarkOrange"), times = 1, len = nrow(chr))

# Filter and organize data for plotting
LoessOut <- LoessOut %>% as_tibble() %>% filter(V1 %in% chr$V1, !is.na(V3))
LoessOut$V1 <- factor(LoessOut$V1, levels = chr$V1)

# Plotting
pbulk <- ggplot(data = LoessOut, aes(x = V2, y = V3)) +
  geom_point(aes(color = V1), size = 0.5) +
  geom_line(aes(x = V2, y = V4), color = "black", linewidth = 1) +
  geom_hline(yintercept = cutvalue, color = "red", linetype = 2, linewidth = 1) +
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

# Save plots
ggsave(paste(outpre, "plot-EDloessAICc.pdf", sep = "."), plot = pbulk, device = "pdf", width = 10,height = 5 )
ggsave(paste(outpre, "plot-EDloessAICc.png", sep = "."), plot = pbulk, device = "png", width = 10,height = 5 )
