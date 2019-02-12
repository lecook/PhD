## Laura E Cook
## Oct 2018
## This script takes the output from MEME FIMO and counts the number of TF matches in the DNA sequences provided to MEME (TWAR sequences)
## It then plots the matches for each TWAR (DNA sequence)
## WORK IN PROGRESS

# Load in packages
library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales) # for ggplot
library(gridExtra) # for ggplot 
library(dplyr)

setwd("~/phd/OneDrive - The University of Melbourne/tfbs_project/results/meme")

options(width=150)

# Read in tables
tcyn <- read.table("tcyn_CFmotifs_HOCOMOCO_core_80%_overlap.tsv", stringsAsFactors=F, sep="\t", header=T)
clup <- read.table("clup_CFmotifs_HOCOMOCO_core_80%_overlap.tsv", stringsAsFactors=F, sep="\t", header=T)
ailMel1 <- read.table("ailMel1_CFmotifs_HOCOMOCO_core_80%_overlap.tsv", stringsAsFactors=F, sep="\t", header=T)
sarHar1 <- read.table("sarHar1_CFmotifs_HOCOMOCO_core_80%_overlap.tsv", stringsAsFactors=F, sep="\t", header=T)

# Clean up motif names
tcyn$motif_id <- gsub(pattern="_H.*", replacement="", x=tcyn$motif_id)
clup$motif_id <- gsub(pattern="_H.*", replacement="", x=clup$motif_id)
ailMel1$motif_id <- gsub(pattern="_H.*", replacement="", x=ailMel1$motif_id)
sarHar1$motif_id <- gsub(pattern="_H.*", replacement="", x=sarHar1$motif_id)

# Clean up sequence names
tcyn$sequence_name <- gsub(pattern="Tcyn_", replacement="", x=tcyn$sequence_name)
clup$sequence_name <- gsub(pattern="Clup_", replacement="", x=clup$sequence_name)
ailMel1$sequence_name <- gsub(pattern="ailMel1_", replacement="", x=ailMel1$sequence_name)
sarHar1$sequence_name <- gsub(pattern="sarHar1_", replacement="", x=sarHar1$sequence_name)

# Add species to datatables
tcyn["species"] <- "tcyn"
ailMel1["species"] <- "ailMel1"
sarHar1["species"] <- "sarHar1"
clup["species"] <- "clup"


# Combine tables, filter by p-value and count number of motifs per species
matches <- rbind(tcyn, clup, ailMel1, sarHar1)
matches$p.value <- format(matches$p.value, scientific=FALSE)
matches_0.01 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.01") %>% count(motif_id, species)
matches_0.001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.001") %>% count(motif_id, species)
matches_0.0001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.0001") %>% count(motif_id, species)
matches_0.00001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.00001") %>% count(motif_id, species)
matches_0.000001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.000001") %>% count(motif_id, species)


# Make plot of binding events for each PWM
pdf("matches_per_PWM.pdf", height=6)
ggplot(matches_0.01, aes(motif_id, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.01	total matches = 142133") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(matches_0.001, aes(motif_id, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.001		total matches = 18520") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(matches_0.0001, aes(motif_id, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.0001	total matches = 2562") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(matches_0.00001, aes(motif_id, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.00001		total matches = 321") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(matches_0.000001, aes(motif_id, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.000001	total matches = 25") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Combine tables, filter by p-value and count number of matches per TWAR
matches <- rbind(tcyn, clup, ailMel1, sarHar1)
matches$p.value <- format(matches$p.value, scientific=FALSE)
matches_0.01 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.01") %>% count(sequence_name, species)
matches_0.001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.001") %>% count(sequence_name, species)
matches_0.0001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.0001") %>% count(sequence_name, species)
matches_0.00001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.00001") %>% count(sequence_name, species)
matches_0.000001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.000001") %>% count(sequence_name, species)


pdf("matches_per_TWAR.pdf", height=6, width=13)
ggplot(matches_0.01, aes(sequence_name, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.01 	~725 matches per TWAR") c
ggplot(matches_0.001, aes(sequence_name, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.001	~95 matches per TWAR") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
ggplot(matches_0.0001, aes(sequence_name, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.0001	~13 matches per TWAR") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
ggplot(matches_0.00001, aes(sequence_name, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.00001	~3 matches per TWAR") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
ggplot(matches_0.000001, aes(sequence_name, n, colour = species)) + geom_point() + ggtitle("matches per PWM < 0.000001	~1.6 matches per TWAR") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
dev.off()




matches <- rbind(tcyn, clup, ailMel1, sarHar1)
matches$p.value <- format(matches$p.value, scientific=FALSE)
matches_0.01 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value <= "0.01") 
matches_0.001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.001")
matches_0.0001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.0001") 
matches_0.00001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.00001") 
matches_0.000001 <- matches %>% select(motif_id, p.value, species, sequence_name) %>% filter(p.value < "0.000001")


# Count number of binding events for each motif for each species
#tcyn_motif_counts <- count(tcyn$motif_id)
#ailMel1_motif_counts <- count(ailMel1$motif_id)
#sarHar1_motif_counts <- count(sarHar1$motif_id)
#clup_motif_counts <- count(clup$motif_id)

# Count number of binding events for each TWAR for each species
#tcyn_TWAR_counts <- count(tcyn$sequence_name)
#ailMel1_TWAR_counts <- count(ailMel1$sequence_name)
#sarHar1_TWAR_counts <- count(sarHar1$sequence_name)
#clup_TWAR_counts <- count(clup$sequence_name)

# Combine tables
#motif_counts <- rbind(tcyn_motif_counts, clup_motif_counts, ailMel1_motif_counts, sarHar1_motif_counts)
#twar_counts <- rbind(tcyn_TWAR_counts, clup_TWAR_counts, ailMel1_TWAR_counts, sarHar1_TWAR_counts)


# Make plot of motif counts
pdf("motif_counts.pdf", height=6)
ggplot(motif_counts, aes(x, freq, colour = species)) + geom_point() + ggtitle("motif_counts") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Subset TWAR count data by chromosome

chr1 <- twar_counts[grep("chr1_.*", twar_counts$x),]
chr2 <- twar_counts[grep("chr2_.*", twar_counts$x),]
chr3 <- twar_counts[grep("chr3_.*", twar_counts$x),]
chr4 <- twar_counts[grep("chr4_.*", twar_counts$x),]
chr5 <- twar_counts[grep("chr5_.*", twar_counts$x),]
chr6 <- twar_counts[grep("chr6_.*", twar_counts$x),]
chr7 <- twar_counts[grep("chr7_.*", twar_counts$x),]
chr8 <- twar_counts[grep("chr8_.*", twar_counts$x),]
chr9 <- twar_counts[grep("chr9_.*", twar_counts$x),]
chr10 <- twar_counts[grep("chr10_.*", twar_counts$x),]
chr11 <- twar_counts[grep("chr11_.*", twar_counts$x),]
chr12 <- twar_counts[grep("chr12_.*", twar_counts$x),]
chr13 <- twar_counts[grep("chr13_.*", twar_counts$x),]
chr14 <- twar_counts[grep("chr14_.*", twar_counts$x),]
chr15 <- twar_counts[grep("chr15_.*", twar_counts$x),]
chr16 <- twar_counts[grep("chr16_.*", twar_counts$x),]
chr17 <- twar_counts[grep("chr17_.*", twar_counts$x),]
chr18 <- twar_counts[grep("chr18_.*", twar_counts$x),]
chr19 <- twar_counts[grep("chr19_.*", twar_counts$x),]

# Save twar count plots into a PDF
pdf("data/matches_per_twar.pdf", width=18, height=16)
p1 <- ggplot(chr1, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr1") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p2 <- ggplot(chr2, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr2") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p3 <- ggplot(chr3, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr3") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p4 <- ggplot(chr4, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr4") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p5 <- ggplot(chr5, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr5") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p6 <- ggplot(chr6, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr6") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p7 <- ggplot(chr7, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr7") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p8 <- ggplot(chr8, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr8") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p9 <- ggplot(chr9, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr9") + xlab("") + ylab("l") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)
p10 <- ggplot(chr10, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr10") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p11 <- ggplot(chr11, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr11") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p12 <- ggplot(chr12, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr12") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p13 <- ggplot(chr13, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr13") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p14 <- ggplot(chr14, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr14") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p15 <- ggplot(chr15, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr15") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p16 <- ggplot(chr16, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr16") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p17 <- ggplot(chr17, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr17") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
p18 <- ggplot(chr18, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr18") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1)) %+replace% theme(legend.position="none")
ggarrange(p10, p11, p12, p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 3)
p19 <- ggplot(chr19, aes(x, freq, colour = species)) + geom_point() + ggtitle("chr19") + xlab("") + ylab("") + theme_bw() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggarrange(p19, ncol = 3, nrow = 3)
dev.off()

# Save motif_count files

