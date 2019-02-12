## Laura E Cook
## 17 Oct 2018
## This script plots the length of DNA sequences from a fasta file

library("ggplot2")
library("dplyr")
library("tidyverse")
library("ggpubr")

# Read in tables

ailMel1 <- read.table("data/fasta_TWAR_files/ailMel1_TWAR_length.txt")
clup <- read.table("data/fasta_TWAR_files/Clup_TWAR_length.txt")
tcyn <- read.table("data/fasta_TWAR_files/Tcyn_TWAR_length.txt")
sarHar1 <- read.table("data/fasta_TWAR_files/sarHar1_TWAR_length.txt")


# Add species column to the dataframe
tcyn["species"] <- "tcyn"
clup["species"] <- "clup"
ailMel1["species"] <- "ailMel1"
sarHar1["species"] <- "sarHar1"

# Clean up TWAR names
tcyn$V1 <- gsub(pattern="Tcyn_", replacement="", x=tcyn$V1)
clup$V1 <- gsub(pattern="Clup_", replacement="", x=clup$V1)
ailMel1$V1 <- gsub(pattern="ailMel1_", replacement="", x=ailMel1$V1)
sarHar1$V1 <- gsub(pattern="sarHar1_", replacement="", x=sarHar1$V1)

# Combine tables
data <- rbind(tcyn, clup, ailMel1, sarHar1)

# Subset data by chromosome

chr1 <- data[grep("chr1_.*", data$V1),]
chr2 <- data[grep("chr2_.*", data$V1),]
chr3 <- data[grep("chr3_.*", data$V1),]
chr4 <- data[grep("chr4_.*", data$V1),]
chr5 <- data[grep("chr5_.*", data$V1),]
chr6 <- data[grep("chr6_.*", data$V1),]
chr7 <- data[grep("chr7_.*", data$V1),]
chr8 <- data[grep("chr8_.*", data$V1),]
chr9 <- data[grep("chr9_.*", data$V1),]
chr10 <- data[grep("chr10_.*", data$V1),]
chr11 <- data[grep("chr11_.*", data$V1),]
chr12 <- data[grep("chr12_.*", data$V1),]
chr13 <- data[grep("chr13_.*", data$V1),]
chr14 <- data[grep("chr14_.*", data$V1),]
chr15 <- data[grep("chr15_.*", data$V1),]
chr16 <- data[grep("chr16_.*", data$V1),]
chr17 <- data[grep("chr17_.*", data$V1),]
chr18 <- data[grep("chr18_.*", data$V1),]
chr19 <- data[grep("chr19_.*", data$V1),]

# Save plots into a PDF
pdf("data/TWAR_length_plots2.pdf", width=18, height=16)
p1 <- ggplot(chr1, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr1") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p2 <- ggplot(chr2, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr2") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p3 <- ggplot(chr3, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr3") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p4 <- ggplot(chr4, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr4") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p5 <- ggplot(chr5, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr5") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p6 <- ggplot(chr6, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr6") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p7 <- ggplot(chr7, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr7") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p8 <- ggplot(chr8, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr8") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p9 <- ggplot(chr9, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr9") + xlab("") + ylab("l") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)
p10 <- ggplot(chr10, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr10") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p11 <- ggplot(chr11, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr11") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p12 <- ggplot(chr12, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr12") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p13 <- ggplot(chr13, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr13") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p14 <- ggplot(chr14, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr14") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p15 <- ggplot(chr15, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr15") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p16 <- ggplot(chr16, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr16") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p17 <- ggplot(chr17, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr17") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p18 <- ggplot(chr18, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr18") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
ggarrange(p10, p11, p12, p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 3)
p19 <- ggplot(chr19, aes(V1, V2, colour = species)) + geom_point() + ggtitle("chr19") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggarrange(p19, ncol = 3, nrow = 3)
dev.off()