#!/usr/bin/env Rscript

## 06 Feb 2019
## Laura E Cook, University of Melbourne
## This script processes the data ouput from MEME FIMO into a useable form
## Takes PWM matches per TWAR sequence and compares the intervals between matches within and between species
## Find PWM matches that are unique to each species and unique to combinations of species
## Count the number of unique matches and plot as a volcano plot for each TWAR
## Decide based on those plots which TWARs and PWM matches are interesting to follow up to compare the scores (do this in another script)
## Also find those that are convergently lost in thylacine and wolf or gained in thylacine and wolf but not present in the outgroup species (do this in another script)


###### ------------------------------- SAVE SESSION INFO, LOAD IN PACKAGES AND DATA ----------------------------------- ######

# Load in libraries
library(plyr)
library(data.table) 
library(dplyr)
library(GenomicRanges)
library(stringr)
library(tidyverse)

# Save session info
writeLines(capture.output(sessionInfo()), "sessionInfo_tfbs_analysis.txt")

# Set working directory
setwd("~/phd/OneDrive - The University of Melbourne/tfbs_project/results/meme/")

# Load in table of craniofacial motif matches from MEME FIMO for all species: Thylacine, Wolf, Panda and Devil
m <- read.table("all_CFmotifs_HOCOMOCO_core_80%_overlap_0.001.tsv", header=TRUE, sep="")


###### ------------------------------- FINDING UNIQUE MATCHES BASED ON INTERVALS ----------------------------------- ######

# Create an empty list for twars
twar.list <- list()

# Loop over every unique twar ID 
for(twar in unique(m$sequence_name)){
print(twar)
# Generate a list with the twar IDs
list.name2 <- as.character(twar)

	# Extract the data for that particular TWAR
	twar_subset <- subset(m, sequence_name==twar)

	# Create an empty list for the TF ids
	tf.list <- list()

	# Loop over every unique TF id within the list of matches for this TWAR
	for(tf in unique(twar_subset$motif_id)){

	# Append the list with all the TF names
	list.name <- as.character(tf)

		# Convert dataframe into a GRanges object to find overlaps and unique sites to compare
		GR <- makeGRangesFromDataFrame(subset(twar_subset, motif_id==tf), keep.extra.columns=TRUE, seqnames.field="sequence_name", start.field="start", end.field="stop", strand.field="strand")
		
		# Extract the interval data according to the metadata "species"
		ailMel1 <- GR[which(elementMetadata(GR)[, "species"] %in% "ailMel1")]
		tcyn <- GR[which(elementMetadata(GR)[, "species"] %in% "tcyn")]
		sarHar1 <- GR[which(elementMetadata(GR)[, "species"] %in% "sarHar1")]
		clup <- GR[which(elementMetadata(GR)[, "species"] %in% "clup")]
		

		# Extract the unique matches for each species
		# Can't work out a way to do this in one go so have to look at thylacine no devil, thylacine no panda and thylacine no wolf and then overlap those lists
		# subsetByOverlaps() looks for the overlapping sequences between two GRanges objects, invert=TRUE means that it does the opposite and outputs the sequence found in thylacine but not in devil
		# Ended up using %over% but this doesn't ignore strand so I might have to go back to subsetByOverlaps
		# Should I exclude strand in determining the overlap? Because The matched sequence can be exactly the same but they don't come up as matched because they're
			# on different strands


		## Unique to thylacine
		thylacine_not_devil <- tcyn[!tcyn %over% sarHar1]

		### In thylacine, not in panda
		thylacine_not_panda <- tcyn[!tcyn %over% ailMel1]

		### In thylacine, not in wolf
		thylacine_not_wolf <- tcyn[!tcyn %over% clup]

		# Overlap the list above to find those matches unique to thylacine
		uT <- as(Reduce(subsetByOverlaps, list(thylacine_not_wolf, thylacine_not_panda, thylacine_not_devil)), "data.frame")

		# If there is data present (there isn't always a unique match for every TF) add a column to the data.frame so I can keep track of which match is unique in which species
		# If the number of rows in the dataframe is greater the zero add the column with species name
		# Repeat this for the other species - perhaps I can do this in a loop if I had a list of the species names?? ASK IRENE
		if(nrow(uT)>0){
			uT$unique <- "thylacine"
		}	

		## Unique to wolf
		wolf_not_devil <- clup[!clup %over% sarHar1]
		wolf_not_panda <- clup[!clup %over% ailMel1]
		wolf_not_thylacine <- clup[!clup %over% tcyn]
		uW <- as(Reduce(subsetByOverlaps, list(wolf_not_thylacine, wolf_not_panda, wolf_not_devil)), "data.frame")
		if(nrow(uW)>0){
			uW$unique <- "wolf"
		}	

		## Unique to panda
		panda_not_devil <- ailMel1[!ailMel1 %over% sarHar1]
		panda_not_thylacine <- ailMel1[!ailMel1 %over% tcyn]
		panda_not_wolf <- ailMel1[!ailMel1 %over% clup]
		uP <- as(Reduce(subsetByOverlaps, list(panda_not_wolf, panda_not_thylacine, panda_not_devil)), "data.frame")
		if(nrow(uP)>0){
			uP$unique <- "panda"
		}	

		## Unique to devil
		devil_not_thylacine <- sarHar1[!sarHar1 %over% tcyn]
		devil_not_panda <- sarHar1[!sarHar1 %over% ailMel1]
		devil_not_wolf <- sarHar1[!sarHar1 %over% clup]
		uD <- as(Reduce(subsetByOverlaps, list(devil_not_wolf, devil_not_panda, devil_not_thylacine)), "data.frame")
		if(nrow(uD)>0){
			uD$unique <- "devil"
		}

		# Look at some other combinations

		## Unique to thylacine and wolf
		uTW <-as(Reduce(subsetByOverlaps, list(thylacine_not_panda, thylacine_not_devil, clup)), "data.frame")
		if(nrow(uTW)>0){
			uTW$unique <- "thylacine_wolf"
		}

		## Unique to thylacine and panda
		uTP <-as(Reduce(subsetByOverlaps, list(thylacine_not_wolf, thylacine_not_devil, ailMel1)), "data.frame")
		if(nrow(uTP)>0){
			uTP$unique <- "thylacine_panda"
		}

		## Unique to thylacine and devil
		uTD <-as(Reduce(subsetByOverlaps, list(sarHar1, thylacine_not_panda, thylacine_not_wolf)), "data.frame")
		if(nrow(uTD)>0){
			uTD$unique <- "thylacine_devil"
		}

		## Unique to wolf and panda
		uWP <-as(Reduce(subsetByOverlaps, list(wolf_not_thylacine, wolf_not_devil, ailMel1)), "data.frame")
		if(nrow(uWP)>0){
			uWP$unique <- "wolf_panda"
		}

		## Unique to wolf and panda
		uWD <-as(Reduce(subsetByOverlaps, list(sarHar1, wolf_not_thylacine, wolf_not_panda)), "data.frame")
		if(nrow(uWD)>0){
			uWD$unique <- "wolf_devil"
		}

		## Unique to panda and devil
		uPD <-as(Reduce(subsetByOverlaps, list(sarHar1, panda_not_thylacine, panda_not_wolf)), "data.frame")
		if(nrow(uPD)>0){
			uPD$unique <- "panda_devil"
		}

		## Unique to thylacine, wolf and devil
		uTWD <-as(Reduce(subsetByOverlaps, list(thylacine_not_panda, sarHar1, clup)), "data.frame")
		if(nrow(uTWD)>0){
			uTWD$unique <- "thylacine_wolf_devil"
		}

		## Unique to thylacine, wolf and panda
		uTWP <-as(Reduce(subsetByOverlaps, list(thylacine_not_devil, ailMel1, clup)), "data.frame")
		if(nrow(uTWP)>0){
			uTWP$unique <- "thylacine_wolf_panda"
		}

		## Unique to wolf, panda and devil
		uWPD <-as(Reduce(subsetByOverlaps, list(wolf_not_thylacine, ailMel1, sarHar1)), "data.frame")
		if(nrow(uWPD)>0){
			uWPD$unique <- "wolf_panda_devil"
		}

		## Unique to thylacine, devil and panda
		uTDP <-as(Reduce(subsetByOverlaps, list(thylacine_not_wolf, ailMel1, sarHar1)), "data.frame")
		if(nrow(uTDP)>0){
			uTDP$unique <- "thylacine_devil_panda"
		}

		## Match is present in all four species
		uTWPD <- as(Reduce(subsetByOverlaps, list(tcyn, ailMel1, sarHar1, clup)), "data.frame")
		if(nrow(uTWPD)>0){
			uTWPD$unique <- "thylacine_wolf_panda_devil"
		}

		# For each TF bind all the combinations into 1 dataframe
		df_temp <- rbind(uT, uW, uD, uP, uTW, uTP, uTD, uWP, uWD, uPD, uTWD, uTWP, uWPD, uTDP, uTWPD)
		
		# Add the TF that the loop is on to the dataframe so I know which TF this data is for
		if(nrow(df_temp)>0){
			df_temp$TF <- tf
		}

			# For each TF in my list add the information from the dataframe (dt_temp) 
			tf.list[[list.name]] <- df_temp
	}

	# Combine all the data for every TF in the list
	final.data <- do.call(rbind, tf.list)

twar.list[[list.name2]] <- final.data

} 

data <- do.call(rbind, twar.list)
write.table(data, file="unique_matches.txt")

###### ------------------------------- COUNTING MATCHES ----------------------------------- ######

# Loop over every unique twar ID 
for(twar in unique(data$seqnames)){

	# Extract the data for that particular TWAR from both the unique matches and the total matches
	twar_subset <- subset(data, seqnames==twar)
	twar_subset_nonU <- subset(m, sequence_name==twar)

	# Create an empty list for T ids
	tf.list <- list()

	# Loop over every unique TF id within the list of matches for this TWAR
	for(tf in unique(twar_subset$TF)){
	
	# Generate a list of the tf names and a list of the unique species
	list.name <- as.character(tf)
					
		# Subset the data by tf for both unique matches and total matches
		tf_subset <- subset(twar_subset, TF == tf)
		tf_twar_subset <- subset(twar_subset_nonU, motif_id == tf)

		# Count occurances in unique matches dataset (final.data)
		uT_count <- sum(str_count(tf_subset$unique, "\\bthylacine\\b"))
		uD_count <-  sum(str_count(tf_subset$unique, "\\bdevil\\b"))
		uP_count <- sum(str_count(tf_subset$unique, "\\bpanda\\b"))
		uW_count <- sum(str_count(tf_subset$unique, "\\bwolf\\b"))
		uTW_count <- sum(str_count(tf_subset$unique, "\\bthylacine_wolf\\b"))
		uTP_count <- sum(str_count(tf_subset$unique, "\\bthylacine_panda\\b"))
		uTD_count <- sum(str_count(tf_subset$unique, "\\bthylacine_devil\\b"))
		uWD_count <- sum(str_count(tf_subset$unique, "\\bwolf_devil\\b"))
		uPD_count <- sum(str_count(tf_subset$unique, "\\bpanda_devil\\b"))
		uWP_count <- sum(str_count(tf_subset$unique, "\\bwolf_panda\\b"))
		uTWD_count <- sum(str_count(tf_subset$unique, "\\bthylacine_wolf_devil\\b"))
		uTWP_count <- sum(str_count(tf_subset$unique, "\\bthylacine_wolf_panda\\b"))
		uWPD_count <- sum(str_count(tf_subset$unique, "\\bwolf_panda_devil\\b"))
		uTDP_count <- sum(str_count(tf_subset$unique, "\\bthylacine_devil_panda\\b"))
		uTWPD_count <- sum(str_count(tf_subset$unique, "\\bthylacine_wolf_panda_devil\\b"))

		# Count all occurances from original dataset (total matches)
		T <- sum(str_count(tf_twar_subset$species, "\\btcyn\\b"))
		W <- sum(str_count(tf_twar_subset$species, "clup"))
		P <- sum(str_count(tf_twar_subset$species, "ailMel1"))
		D <- sum(str_count(tf_twar_subset$species, "sarHar1"))


		# Save these into a new df
		df_count <- data.frame(T, W, P, D, uT_count, uW_count, uD_count, uP_count, uTW_count, uTP_count, uTD_count, uWD_count, uPD_count, uWP_count, uTWD_count, uTWP_count, uWPD_count, uTDP_count, uTWPD_count)
			# Sort by list of tfs
			tf.list[[list.name]] <- df_count

		# Combine all the count data for every TF in the list
		count.data <- do.call(rbind, tf.list)	
	
	}
	# Ouput a new table for every TWAR 
	write.table(count.data, paste(twar, ".txt", sep=""))


} 








