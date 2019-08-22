#-----------------------------------------------------------------------------------------------------------------------
# Program description
# Program to align oligonucleotide sequences by conserved sequence motifs and group into related families.
# Program runs fasta files 
# Program outputs a txt file with path and name input by user

#------------------------------------------------------------------------------------------------------------------------

# User input values. 
# Enter values on right-hand side of arrows below, press command A, and click run above to analyze sequences.

# Enter filename. Use * as a wild card to process multiple files in a single run
# fileNames <- "NP*.fasta"
fileNames <- "NP01-3-A01_R1_001.fastq_merged.fasta"

# Enter path to files 
path <- "/Users/VelocityCompute/Documents/Plate3/merged"
# path <- "/Users/VelocityCompute/Documents"

# Input path for alignment report 
report_path <- "/Users/VelocityCompute/Documents/Plate3/alignments"
# report_path <- "/Users/VelocityCompute/Documents"

# Input desired length of conserved motif
motif.length <- 15

# Input allowed number of insertions, deletions, and substitutions (0 = perfect match)
string_dif_threshold <- 1

# Input minimum sequence family size desired
family_size <- 3

# Input fixed regions to be deleted from each sequence
# Input 5pai for the sense strands, the longest sequence on the 5' end of library A (internal primer)
fixed1 <- 'AATACGACTCACTATAGGGAGAGGCGAAAAGAGCCACAAAGCA'

# Input length of fixed1
fixed1_length <- 44

# Input the  3' primer for the sense strands
# fixed2 <- 'GAAGGAGGACGAGGAAGAACACGACAGGC'
fixed2 <- 'GAAGGAGGACGAGGAAGAACACGAC'
          
# Input length of fixed2
fixed2_length <- 29

# Input the 5' primer for the antisense strands (reverse complement of fixed1)
fixed3 <- 'TGCTTTGTGGCTCTTTTCGCCTCTCCCTATAGTGAGTCGTATTA'

# Include bead binders with alignment (enter TRUE or FALSE  in all caps)
# This file is a list of non-replicated bead binders that is co-aligned with the non-replicated target sequences
include_bead_binders <- TRUE
# Enter path and filename for bead binders
bead_background <- "/Users/VelocityCompute/Documents/Bead Files/NP08-9-8_R1_001_trimmed_BeadBinders.txt"

# Enter path and filename for the complete bead background sequence file
# This is the complete file (replicates, antisense and sense) of bead binders that is used to check for matches of the winning target sequences
complete_bead_background_file <- "/Users/VelocityCompute/Documents/Bead Files/NP08-9-8_R1_001_beaddb_merged.fasta"

# Enter true to append summary statistics to csv file
stat_append <- FALSE
# Enter path and file name for stats
stat_file <- "/Users/VelocityCompute/Documents/Plate3/AlignmentStats-Plate_3.csv"

# Enter true to append sequences to csv file
append_picks <- FALSE
# Enter a file to export top 10 most replicated sequences that appear in alignment
pick_file <- '/Users/VelocityCompute/Documents/Plate3/Plate 3-picks_from_alignments.csv'

# To limit the number of sequences analyzed, enter TRUE and the number of sequences to analyze
# Most useful for doing a quick alignment, particularly on early round data
# Entire sequence file will be processed to remove duplicates, shortmers, and longmers, but not close matches
# Then number_to_analyze will be selected and close matches will be removed and sequences aligned
# Close matches will not be collected as a data point in the summary statistics output
# Note if number of sequences after filtering > 3000, dataframe will automatically be reduced to the number in number_to_analyze (see line 336)
select_subset <- FALSE
number_to_analyze <- 3000

# Export list of non-bead binders as text file in fasta format
# Enter TRUE in all caps to export list
nonbeadlist <- FALSE

# Export filtered data to fasta text file for use in MEME
# Enter TRUE or FALSE in all caps 
export_fasta <- FALSE


# End of user input section
#------------------------------------------------------------------------------------------------------------------------

# Install and load required packages

# source("http://www.bioconductor.org/biocLite.R")
# install.packages("iterators")
# library(msa)
# system.file("tex", "texshade.sty", package = "msa")
# Loading fasta files and creating a dataframe. (File may contain multiple sequences.) (see: https://stackoverflow.com/questions/21263636/how-to-read-fasta-into-dataframe-and-extract-subsequences-of-fasta-file-in-r)
# library("Biostrings")
# library("iterators")
# install.packages('plyr')
# library(plyr)
# library('DECIPHER')
# install.packages('rlist')
# library(rlist)
# library(stringr)
# install.packages('dplyr')
# library('dplyr')
# install.packages('seqinr')
# library('seqinr')

#-----------------------------------------------------------------------------------------------------
# Code. No user input required below!

# Read in all files 
# Join file_names and path
path_and_file <- file.path(path, fileNames)

# Put all files into list to run alignment on multiple files in loop (glob files must be in list to loop)
file <- as.list(Sys.glob(path_and_file))
# print(file)

# Create loop through all files that glob placed in file list
for(j in file){
  # Read in fasta file to be analyzed
  fastaFile <- readDNAStringSet(j)
  
  # Create path and file for alignment report
  split_fileName <- basename(j)
  alignment_fileName <- paste0("Aligned-MotifLen=", motif.length, "INDELS=", string_dif_threshold, split_fileName)
 
  # Join path and file name for alignment report for each file to be analyzed (used in the open sink command below)
  report_path_and_file <- file.path(report_path, alignment_fileName)

  options(stringsAsFactors = FALSE)

  # Assign variable to sequence name
  seq_name = names(fastaFile)

  # Assign variable to sequence
  sequence = paste(fastaFile)

  # Create dataframe from sequence name and sequence data in the fasta file
  sequence.df.initial <- data.frame(seq_name, sequence)

  # Get total # of sequences in dataframe
  allSeqs <- nrow(sequence.df.initial)

  # Start a timer
  start_time <- Sys.time()

  # Remove rows with "N" in the sequence column
  # -grep acts weird on some dataframes (it deleted all rows on a test fasta data set once)
  # sequence.df <- sequence.df.initial[-grep("N", sequence.df.initial$sequence),]
  sequence.df <- sequence.df.initial[!grepl("N", sequence.df.initial$sequence),]
  # I have now seen !grepl delete all rows as well so may need to use 
  # sequence.df <- sequence.df[grep("N", sequence.df.initial$sequence, invert = TRUE), ]

  # Calculate number of sequences containing Ns
  N_in_seq <- nrow(sequence.df)
  number_of_N <- allSeqs - N_in_seq
  
  # Remove short sequences that can get through if they are around 40 bases long and don't have correct primers to be filtered out
  # Don't remove sequences longer than the expected 113 here though because the 40N random region can be extracted from them
  sequence.df <- sequence.df[!nchar(as.character(sequence.df$sequence)) < 108, ]
  
  # Calculate number of sequences in sequence.df less than 108
  first_shortmer_cut <- N_in_seq - nrow(sequence.df)

  # Convert antisense strands to sense strands by converting them to their reverse complements
  # Grep will only look for exact matches of the fixed region so sequences with indelsubs in the fixed region are lost
  # Partition antisense strands into a new df by searching for fixed3
  # Subset df by fixed3 to keep the antisense strands and remove the sense strands
  # Grep will subset by exact matches, agrep by the input number of mismatches
  # sequence.df.antisense <- sequence.df[grep(fixed3, sequence.df$sequence),]
  sequence.df.antisense <- sequence.df[agrep(fixed3, sequence.df$sequence, max.distance = 2, value = FALSE ),]

  # Get # of antisense before removing duplicates
  anti_before <- nrow(sequence.df.antisense)

  # Remove exact duplicates
  sequence.df.antisense.nodups <- sequence.df.antisense[!duplicated(sequence.df.antisense[2]), ]

  # Get # of unique (no exact matches) antisense strands
  no_of_antisense <- nrow(sequence.df.antisense.nodups)

  # Calculate # of antisense duplicates
  no_of_antisenseDups <- anti_before - no_of_antisense

  # Convert antisense sequences to their reverse complements
  sequence.df.antisense.nodups$comp <- sapply(sequence.df.antisense.nodups$sequence, function(x) as.character(reverseComplement(DNAString(x))))

  # Delete original sequences
  sequence.df.antisense.nodups <- sequence.df.antisense.nodups[, -2]

  # Rename antisense sequence column to match sense dataframe so antisense and sense dataframes can be merged
  colnames(sequence.df.antisense.nodups)[2] <- "sequence"

  # Subset df by fixed1 to keep the sense strands and remove the antisense strands
  # Use grep for exact matches to fixed1 and agrep for mismatches in fixed1
  # sequence.df.sense <- sequence.df[grep(fixed1, sequence.df$sequence),]
  sequence.df.sense <- sequence.df[agrep(fixed1, sequence.df$sequence, max.distance = 2, value = FALSE ),]

  # Get # of sense strands before subtracting duplicates
  sense_before <- nrow(sequence.df.sense)

  # Remove exact duplicates
  sequence.df.sense.nodups <- sequence.df.sense[!duplicated(sequence.df.sense[2]), ]

  # Get # of sense strands after subtracting out exact matches
  no_of_sense <- nrow(sequence.df.sense.nodups)

  # Calculate number of sense duplicates
  no_of_senseDups <- (nrow(sequence.df.sense)) - no_of_sense

  # Merge sense and antisense dataframes into a dataframe of sequences containing no perfect matches
  sequence.df.noperfectmatches.merged <- bind_rows(sequence.df.sense.nodups, sequence.df.antisense.nodups)
  sense_antisense_count <- (nrow(sequence.df.noperfectmatches.merged))

  # Remove exact duplicates again, just in case there were duplicates across sense and antisense sets
  sequence.df.noperfectmatches <- sequence.df.noperfectmatches.merged[!duplicated(sequence.df.noperfectmatches.merged[2]), ]
  seq_number <- nrow(sequence.df.noperfectmatches)
  moreDups <- sense_antisense_count - seq_number

  # Calculate total number of duplicate sequences sense dups + antisense dups
  # This does count sequences in the sense df that were also in the antisense df
  no_of_dups <- no_of_senseDups + no_of_antisenseDups + moreDups
  
  #-------------------------------------------------------------------------------------------------------------------------
  # Several strategies for removing the 5' and 3' primer regions are presented below
  
  # Remove 5' primer and everything that comes before and 3' primer and everything that comes after
  # Caveat to this approach is that primers with snps won't be found and sequence will be excluded
  # Remove everything before and including the 5' primer
  # Leave first base off of fixed1 because it is frequently absent in the sequencing data
  # sequence.df.noperfectmatches$sequence <- gsub(".*AATACGACTCACTATAGGGAGAGGCGAAAAGAGCCACAAAGCA", "", sequence.df.noperfectmatches$sequence)
  # End with GTA as a number of 5' primers end with it instead of GCA
  # sequence.df.noperfectmatches$sequence <- gsub(".*AATACGACTCACTATAGGGAGAGGCGAAAAGAGCCACAAAGTA", "", sequence.df.noperfectmatches$sequence)
  # Remove 3' primer and everything after (note last several bases have been ommitted because they often do not appear in the sequences)
  # sequence.df.noperfectmatches$sequence <- gsub("GAAGGAGGACGAGGAAGAACA.*", "", sequence.df.noperfectmatches$sequence)
  
  # Shorten the 5' primer sequence to remove primer and all bases upstream
  # Make the last CA wildcards because they seem to vary quite a bit
  sequence.df.noperfectmatches$sequence <- gsub(".*AGGCGAAAAGAGCCACAAAG..", "", sequence.df.noperfectmatches$sequence)
  # Shorten the 3' primer sequence to remove primer and bases downstream  
  sequence.df.noperfectmatches$sequence <- gsub("GAAGGAGGACGAGGAAGAACA.*", "", sequence.df.noperfectmatches$sequence)

  # Alternate strategy: Remove 1st 44 bases (length of sense 5' primer) and last 29 bases(length of sense 3' primer)
  # Best approach unless long tails appear after 3' fixed region on the sequences
  # Not good for primer dimers
  # Make sure you know the lengths of the primers!
  # Potential problem is that the fixed regions might not be in the same place in every sequence
  # Remove first 44
  # str_sub(sequence.df.noperfectmatches$sequence, 1,fixed1_length) <- ""
  # Remove last 29
  # str_sub(sequence.df.noperfectmatches$sequence, -fixed2_length) <- ""
#---------------------------------------------------------------------------------------------------------------------------
  # Check sequence length
  # length_check <- sequence.df.noperfectmatches[nchar(as.character(sequence.df.noperfectmatches$sequence)) > 44, ]

  # Sometimes very short sequences end up in sequence.df.nodups (even 2mers!).
  # This is likely the result of sequences in the original database that contain fixed regions + a short-mer.
  # Thus, delete sequences shorter than 37 here when looking for 40mers or 47 when looking for 50mers (e.g., random region + 5 primer bases on each   side of random of random region)
  sequence.df.noperfectmatches <- sequence.df.noperfectmatches[!nchar(as.character(sequence.df.noperfectmatches$sequence)) < 37, ]

  # Get number of shortmers, which is the sum of the number shorter than 108 above (before removing primers) and the number <37 or 47
  current_total <- nrow(sequence.df.noperfectmatches)
  second_shortmer_cut <- seq_number - current_total
  shortmers <- second_shortmer_cut + first_shortmer_cut

  # Delete sequences longer than 43
  sequence.df.noperfectmatches <- sequence.df.noperfectmatches[!nchar(as.character(sequence.df.noperfectmatches$sequence)) > 43, ]

  # Get number of sequences longer than 43
  longmers <- current_total - nrow(sequence.df.noperfectmatches)

  # Sort sequence.df.noperfectmatches from longest sequence to shortest
  # This is helpful in filtering out close matches using while loop below
  # Without this, some long sequences get deleted in favor of shorter sequences
  sequence.df.noperfectmatches <- sequence.df.noperfectmatches[order(nchar(sequence.df.noperfectmatches$sequence), decreasing = TRUE), ]

  # Get current # of sequences in sequence.df.noperfectmatches
  no_after_dups <- nrow(sequence.df.noperfectmatches)

# If select_subset = TRUE, and the number of sequences to analyze entered < the number currently in sequence.df.nodups, then execute
# This will save time on the deletion of mismatches below
if (select_subset){
  # Subset df into x sequences chosen randomly
  if (nrow(sequence.df.noperfectmatches) > number_to_analyze){
    sequence.df.noperfectmatches <- sequence.df.noperfectmatches[sample(nrow(sequence.df.noperfectmatches), number_to_analyze), ]
    # Or choose the first x sequences of the dataframe (good for comparing runs and testing code)
    # sequence.df.noperfectmatches <- sequence.df.noperfectmatches[1:number_to_analyze, ]
  }
}

# Search through sequences and remove sequences that are similar to within 4 bases
# Create empty dataframe with columns seq_name and sequence
sequence.df.nodups <- data.frame("seq_name" = character(0), "sequence" = character(0))

# Start a counter for the row # in the new dataframe
i = 1
# Create a list for the matching index values in sequence.df.noperfectmatches
matchRowList <- list()
# While loop to check to see if sequence.df.noperfectmatches is empty
# If not empty, take the entry in row 1 and move it to sequence.df.nodups
# Compare the newest entry in sequence.df.nodups to sequence.df.noperfectmatches and delete close matches
# Helps to start with the longest sequence in the dataframe, which is why noperfectmatches was ordered from longest to shortes
# Stop when sequence.df.noperfectmatches is empty
while(nrow(sequence.df.noperfectmatches) > 0){
  # Copy row from sequence.df.noperfectmatches to sequence.df.nodups
  sequence.df.nodups[i, ] <- sequence.df.noperfectmatches[1, ]
  
  # Compare sequence in row 1 of sequences.df.nodups to all sequences in sequence.df.noperfectmatches to find partial matches
  # Returns index values when Value = True
  matchRows <- agrep(sequence.df.nodups[i, 2], sequence.df.noperfectmatches$sequence, max.distance = 4, value = FALSE)
  
  # Put index value in matchRowList
  matchRowList[["index"]] <- matchRows
  # print (matchRowList)
  # Delete all rows from sequence.df.noperfectmatches that appear in the matchRowList
  for (j in matchRowList){
    sequence.df.noperfectmatches <- sequence.df.noperfectmatches[-(j), ]
    }
    i = i + 1
  }


# Calculate # of sequences removed because of up to 4 base mismatches
# If only a fraction of dataset was analyzed then closeMatches won't be calculated 
if (select_subset){
  closeMatches <- "Not Calculated Due to Data Subsetting"
  }else{
  closeMatches <- no_after_dups - (nrow(sequence.df.nodups))
  }

# If the number of sequences that made it through all filters is >3000, take a subset (set by number_to_analyze at top) of the sequences forward to save time
if (nrow(sequence.df.nodups) > number_to_analyze){
  sequence.df.nodups <- sequence.df.nodups[sample(nrow(sequence.df.nodups), number_to_analyze), ]
  # Set a boolean that triggers print statement below indicating subsetting has occurred
  subbed <- TRUE
  }else{
  subbed <- FALSE
  }

# Get total number of sequences to be aligned
no_for_alignment <- nrow(sequence.df.nodups)

# Load bead binder sequences and merge with protein target sequences
if (include_bead_binders){
  # Input path and file name of fasta file
  beadFile <- readDNAStringSet(bead_background)
  # Create dataframe
  # Assign variable to sequence name
  seq_name = names(beadFile)
  # Assign variable to sequence
  sequence = paste(beadFile)
  # Create dataframe from sequence name and sequence data in the fasta file
  bead.df <- data.frame(seq_name, sequence)
  # Merge the bead and protein dataframes (bead.df and sequence.df.nodups)
  sequence.df.nodups <- bind_rows(sequence.df.nodups, bead.df)
}

# Convert sequence elements to strings
sequence.df.nodups$sequence <- as.character(sequence.df.nodups$sequence)

# Export text file filtered sequence file (sequence.df.nodups)
# Open the sink to generate fasta text file from print statements below
# This will write a file that MEME can import, but R will not 
if (export_fasta){
  # Create path and file for fasta
  split_fileName <- basename(j)
  filtered_fileName <- paste0("Filtered-", split_fileName)
  
  # Join path and file name for alignment report for each file to be analyzed (used in the open sink command below)
  filtered_path_and_file <- file.path(path, alignment_fileName)
  
  sink(filtered_path_and_file, type = c("output", "message"))
  for (row in 1:nrow(sequence.df.nodups)){
    cat(paste0(">", sequence.df.nodups[row, 1], "\r"))
    cat(paste0(sequence.df.nodups[row, 2], "\r"))
  }
  # Close the sink
  sink()
}

#--------------------------------------------------------------------------------------------------------------------------
# Use to convert the filtered sequences to a fasta file that R can read
# if (export_fasta){
# Create path and file for fasta
# split_fileName <- basename(j)
# filtered_fileName <- paste0("Filtered-", split_fileName)
# Join path and file name for alignment report for each file to be analyzed (used in the open sink command below)
# filtered_path_and_file <- file.path(path, alignment_fileName)
# nodups.names <- as.list(sequence.df.nodups$seq_name)
# nodups.sequences <- as.list(sequence.df.nodups$sequence)
# write.fasta(sequences = nodups.sequences, names = nodups.names, file.out = filtered_path_and_file, nbchar = 80)
#}
#--------------------------------------------------------------------------------------------------------------------------

# Loop to create motifs for all sequences in the dataframe
# Create empty list for motifs
motifList <- list()

# Loop through all sequences, chopping each sequence into pieces of motif.length
for (i in 1:nrow(sequence.df.nodups)){
  leftlims <- 1:(nchar(sequence.df.nodups[i, 2]) - (motif.length - 1))
  rightlims <- motif.length:nchar(sequence.df.nodups[i, 2])
  substrings <- mapply(substr, sequence.df.nodups[i, 2], leftlims, rightlims, USE.NAMES = FALSE)
  # print(substrings) # Prints all possible n-mer motifs for every sequence in dataframe
  name <- paste('item:', i, sep ='')
  motifList[[name]] <- substrings
}
# print (motifList)

# Determine length of all elements in motifList in order to create a dataframe
numberofRows <- sum(sapply(motifList,length))
# print (numberofRows)

# Create dataframe from motif list
# Motif list is a huge list of lists. Creating a dataframe separates each motif in its own row of the df.
motif.df <- data.frame(matrix(unlist(motifList), nrow = numberofRows, byrow = T), stringsAsFactors = FALSE)

# Create list of unique conserved motifs to remove duplicates
motifs.unique <- unique(motif.df)
# print (motifs.unique)

# Convert list to character vector because agrep requires the pattern to be a character vector
motifs.unique.vector <- unlist(motifs.unique, use.names = FALSE)
# print (motifs.unique.vector)

# Determine families by conserved motif: uses Levenshtein edit distance to allow for indels and mutations.
# max.distance is the total number of allowed substitutions, deletions, and insertions
# Create list for sequence familes.
seqFamilyList <- list()
for (i in motifs.unique.vector){
  seqFamily <- sequence.df.nodups[sapply(i, agrep, sequence.df.nodups$sequence, max.distance = string_dif_threshold), ]
  
# Append sequence family and conserved motif to list.
  name <- paste0(i, sep = '')
  seqFamilyList[[name]] <- seqFamily$sequence
}

# Keep elements of seqFamilyList with > n elements in the family where n is the desired family size
new_seqFamilyList <- (seqFamilyList)[lengths(seqFamilyList) >= family_size]
# print (new_seqFamilyList)

# Remove Duplicate families if they exist. 
# That is, if the same sequences are in two different families, one family will be deleted!
# This will keep the list names (conserved motifs).
nodups_seqFamilyList <- new_seqFamilyList[!duplicated(new_seqFamilyList)]
# print (nodups_seqFamilyList)

#------------------------------------------------------------------------------------------------------

# Add stars to the conserved regions of every sequence 

# Extract conserved motifs from nodups_seqFamilyList
match_patterns <- names(nodups_seqFamilyList)

# Count the patterns as we print them out
pat_cnt <- 0

# Character to use to highlight the sequences found
hilite_chr <- "*"

# Carriage return\Line feed
crlf <- "\r\n"

# Get number of rows
nrows <- nrow(sequence.df.nodups)

# print ("WARNING: DOES NOT Identify MULTIPLE MATCHES PER LINE")

# Print the threshold being used for string matching tolerance
# cat (paste("String matching tolerance [sum of all insertions, deletions, substitutions for match]:", string_dif_threshold, crlf))

tot_match_cnt <- 0 # count across all seq          

# Create an output file for the alignment.
# Sink prints the contents of the console so however the output looks in the console is how it will look in to txt file.
# Thus, the enviroment window at right may need to be narrowed to print everything on one line. 
sink(report_path_and_file, type = c("output", "message"))

# Loop through conserved motifs
for (pat in match_patterns) {
  pat_cnt <- pat_cnt + 1
  fam_match_cnt <- 0 # count of pattern matches within a given family           
  # Print the pattern and pattern count for the user
  cat (paste0("Family (", pat_cnt, "): ", pat, crlf))
  # Loop through rows for now (might be a bit slower but easier to see what is happening)
  for (row in seq(1:nrows)) {
    # Get the current sequence line from the dataframe 
    seq_txt <- sequence.df.nodups$sequence[row]
    seq_nam <- sequence.df.nodups$seq_name[row]
    # Find partial matches
    # CAUTION: NEED TO EXPAND THIS APPROACH TO MAKE IT WORK FOR MULTIPLE MATCHES IN A LINE
    result <- adist(pat, seq_txt, counts = TRUE, partial = TRUE)
    # Get the match counts
    counts <- attr(result,"counts")
    # Get the total string distance, only use if less than some threshold
    count_sum <- sum(counts)
    # Get the match position offsets
    offsets <- attr(result, "offsets")
    firsts <- offsets[,,1]
    lasts <- offsets[,,2]
    # Make a text string to put under the current sequence that has * from first to last positions
    nrep <- as.numeric(lasts) - as.numeric(firsts)
    #print(nrep)
    if (count_sum <= string_dif_threshold) {
      # Check to see how many matches in this row
      tot_match_cnt <- tot_match_cnt + 1
      fam_match_cnt <- fam_match_cnt + 1
      md_str <- paste0(" INDELSUBS: ", count_sum, " SeqNam: ", seq_nam)
      # Uncomment will leave out match distance from print statement below
      #md_str <- ""
      cat (paste0(seq_txt, md_str, crlf))
      hilite_str <- paste0(str_pad(hilite_chr, firsts, pad = " "),
                           paste(rep(hilite_chr, times = nrep), collapse = "" ))
      cat (paste0(hilite_str, crlf))
      
      # Print the string matching distances
      cat (crlf)
      # Get the match string
      match_str <- substr(x = seq_txt, start = firsts, stop = lasts)
      row_df <- data.frame(seq_pat = pat,
                           seq_str = seq_txt,
                           seq_nam = seq_nam,
                           match_dst = count_sum,
                           match_str = match_str,
                           match_first = firsts,
                           match_last = lasts)
      row.names(row_df) <- tot_match_cnt
      # Append the match to the output data frame
      if (tot_match_cnt == 1) {
        out_df <- row_df
      } else {
        out_df <- rbind(out_df, row_df)
      }  
    } else {
      # Greater than the threshold for distance
      #hilite_str <- paste(rep(" ", times=nchar(seq_txt)), collapse = "")
      nomatchrow <- TRUE
    }  
  } # end of family loop
  # fam cnt here
  cat(paste(fam_match_cnt, "family member(s) found", crlf))
  cat(paste(rep("=", times = 140), collapse = ""), crlf, crlf)
} # end of pattern loop

# Print total number of sequences found in all families
# cat (paste0(tot_match_cnt, " Total matching sequence(s) found across all families", crlf))
# cat(paste0("Does not include multiple matches of the the same pattern in a single sequence", crlf))

# Summarize by pattern
# First make sure at least one family of the desired length was found; ends at line 686
if(exists('out_df') && dim(out_df)[1] > 0){
tmp_df <- out_df[, c("seq_pat")]
tbl <- table(tmp_df)
cat(crlf)
cat("Number of sequences/family:")
cat(crlf)
print(tbl)
cat(crlf)

# Print list of all sequences in families ordered by sequence instead of family
# Arrange out_df in terms of an ordered family tree
fam_tree <- out_df

# Order the dataframe by sequence (seq_str)
fam_tree_ordered <- fam_tree[order(fam_tree$seq_str), ]

# Aggregate all columns by seq_pat
all_data <- aggregate(fam_tree_ordered[ , c(1,3)], by = list(fam_tree_ordered$seq_str), paste, collapse = ",")

# Print each sequence and all families to which those sequences belong
cat(crlf, crlf)
# cat("List of Sequences Belonging to Multiple Families by Sequence:")
cat(crlf)
# Print motif and sequence ID
# for (i in 1:nrow(all_data)){
  # cat("Sequence: ")
  # cat(all_data$seq_str[i], collapse = '\n')
  # cat("Sequence ID: ")
  # cat(all_data$seq_nam[i], collapse = '\n')
  # cat("Conserved Motifs: ")
  # cat(all_data$seq_pat[i], collapse = '\n')
  # cat(crlf)
  # }

# Print list of sequences and motifs for the protein binders and not the bead binders
# Furthermore, do not show a sequence in a bead family
cat("Sequence List:")
cat(crlf)
# Subset fam_tree_ordered into bead sequences
fam_tree_ordered_beadseqs <- fam_tree_ordered[grep("Bead", fam_tree_ordered$seq_nam), ]

# Take bead sequences out of fam_tree_ordered
fam_tree_ordered_nobeadseqs <- fam_tree_ordered[grep("Bead", fam_tree_ordered$seq_nam, invert = TRUE), ]

# Compare motifs in fam_tree_ordered (bead motifs) with motifs in fam_tree_ordered_nobeadseqs and remove rows that match
fam_tree_ordered_nobeadseqs <- fam_tree_ordered_nobeadseqs[!(fam_tree_ordered_nobeadseqs$seq_pat %in% fam_tree_ordered_beadseqs$seq_pat),]

# List each sequence and its conserved motif in order, provided the sequence is not also found in a bead family
if (dim(fam_tree_ordered_nobeadseqs)[1] > 0){
rowTotal <- nrow(fam_tree_ordered_nobeadseqs)
for (row in seq(1:rowTotal)) {
  # Get the current sequence line from the dataframe 
  motif <- fam_tree_ordered_nobeadseqs$seq_pat[row]
  full_sequence <- fam_tree_ordered_nobeadseqs$seq_str[row]
  seq_id <- fam_tree_ordered_nobeadseqs$seq_nam[row]
  # Find partial matches between motif and full_sequence
  matches <- adist(motif, full_sequence, counts = TRUE, partial = TRUE)
  
  # Get the match position offsets
  offset <- attr(matches, "offset")
  match_start <- offset[,,1]
  match_end <- offset[,,2]
  # Make a text string to put under the current sequence that has * from first to last positions
  nRepeats <- as.numeric(match_end) - as.numeric(match_start)
  # Print motif, sequence, and seq id
  cat (paste0(motif, crlf, seq_id, crlf, full_sequence, crlf))
  hilite_str <- paste0(str_pad(hilite_chr, match_start, pad = " "),
                       paste(rep(hilite_chr, times = nRepeats), collapse = "" ))
  cat (paste0(hilite_str, crlf))
  cat(crlf)
  }
}else{
  cat(paste("All protein sequences were found in families containing bead sequences.", crlf, crlf))
# End of if (dim(fam_tree_ordered_nobeadseqs) > 0 statement
}

# Print list of sequences that appear in families without bead sequences
# Collapse the fam_tree dataframe by seq_pat (motif)
fam_tree_byMotif <- aggregate(fam_tree[ , 2:3], by = list(fam_tree$seq_pat), paste, collapse = " ") 

# Subset rows from fam_tree_byMotif that contain the word "Bead" in the sequence name
# This contains all families that contain target sequences along with bead sequences (and families with just bead seqs)
fam_tree_beadSub <- fam_tree_byMotif[grep("Bead", fam_tree_byMotif$seq_nam), ]

# Expand the fam_tree_beadSub dataframe back out to give individual sequences in a vector that were in a family with a bead binder
fam_tree_beadflat <- fam_tree_beadSub$seq_str
fam_tree_BB <- unlist(fam_tree_beadflat, recursive = FALSE)
fam_tree_BB <- strsplit(fam_tree_BB, " ")
fam_tree_BB <- unlist(fam_tree_BB, recursive = FALSE)

# Convert to dataframe for use with grep
fam_tree_BB.df <- data.frame(fam_tree_BB, stringsAsFactors = FALSE)

# Collect sequences from fam_tree that are not in the bead binder dataframe
noBB.df <- fam_tree[!(fam_tree$seq_str %in% fam_tree_BB.df[,1]), ]

# Remove duplicate sequences
noBB.df.nodups <- noBB.df[!duplicated(noBB.df[2]), ]

#--------------------------------------------------------------------------------------------------------------
# Chunk to determine the number of exact and close matches of winning sequences that are present in the original 
# target and bead sequencing files

# Check to see if noBB.df was empty
if(dim(noBB.df)[1] > 0){
  
  # Create a new dataframe from the non-replicated sense sequences from noBB.df
  noBeadBinders.df <- as.data.frame(noBB.df.nodups[, 2], stringsAsFactors = FALSE)
  
  # Sort sequences from longest to shortest
  noBeadBinders.df  <- as.data.frame(noBeadBinders.df[order(nchar(noBeadBinders.df[ , 1]), decreasing = TRUE), ])
  
  # Remove closely matching sequences
  # Create new dataframe to place non-closely matching sequences
  noBB_sequences.df <- data.frame("sequence" = character(0))
  
  # Start a counter for the row # in the new dataframe
  a = 1
  # Create a list for the matches as index values
  newBBList <- list()
  # While loop to check to see if noBeadBinders.df  is empty
  # If not empty, take the entry in row 1 and move it to noBB_sequences.df 
  # Compare the newest entry in noBB_sequences.df  to noBeadBinders.df and delete close matches
  # Helps to start with the longest sequence in the dataframe, which is why noBeadBinders.df was ordered from longest to shortest
  # Stop when noBB_sequences.df is empty
  while(nrow(noBeadBinders.df) > 0){
    # Copy row from noBeadBinders.df to noBB_sequences.df
    noBB_sequences.df[a, 1] <- noBeadBinders.df[1, 1]
    
    # Compare sequence in row 1 of noBB_sequences.df to all sequences in noBeadBinders.df to find partial matches
    # Returns index values when Value = True
    match <- agrep(noBB_sequences.df[a, 1], noBeadBinders.df[, 1], max.distance = 3, value = FALSE)
    # Put index value in newBBList
    newBBList[["index"]] <- match
    
    # Delete all rows from noBeadBinders that appear in the newBBList
    for (b in newBBList){
      noBeadBinders.df <- as.data.frame(noBeadBinders.df[-(b), 1])
    }
    a = a + 1
  }
  
  # Change column name
  colnames(noBB_sequences.df )[1] <- "sense_sequence"
  
  # Add a column in which every sense sequence has been converted to antisense
  noBB_sequences.df$antisense_sequence <- sapply(noBB_sequences.df$sense_sequence, function(x) as.character(reverseComplement(DNAString(x))))
  
  # Count the number of sense and anatisense sequences in the original target and bead datasets with matches to within 4 bases
  # For each sense sequence in noBB_sequences.df, count the number of partial matches to each sequence in sequence.df and in complete_bead_background_file
  
  # Read in bead background file
  bead_file <- complete_bead_background_file
  path_and_file <- file.path(bead_file)
  fastaFile <- readDNAStringSet(path_and_file)
  
  options(stringsAsFactors = FALSE)
  
  # Assign variable to sequence name
  seq_name = names(fastaFile)
  
  # Assign variable to sequence
  sequence = paste(fastaFile)
  
  # Create dataframe from sequence name and sequence data in the fasta file
  bead_sequence.df <- data.frame(seq_name, sequence)
  
  # Loop through files looking for matches
  for (row in 1:nrow(noBB_sequences.df)){
    # Start a counter for number of partially matching sense sequences in target sequence dataframe
    i = 0
    # Start a counter for number of partially matching antisense sequences in target sequence dataframe
    j = 0
    # Start a counter for number of partially matching sense sequences in bead sequence dataframe
    p = 0
    # Start a counter for number of partially matching antisense sequences in bead sequence dataframe
    q = 0
    
    for (seq in 1:nrow(sequence.df)){
      # Use agrepl to return a logical vector for a match to the target sequence file and then count number of 'True'
      sense_matchCount <- agrepl(noBB_sequences.df[row, 1], sequence.df[seq, 2], max.distance = 1)
      antisense_matchCount <- agrepl(noBB_sequences.df[row, 2], sequence.df[seq, 2], max.distance = 1)
      
      # Use agrepl to return a logical vector for a match to the bead sequence file and then count number of 'True'
      # Add up the number of True responses for each sequence and append to sense_countList
      sense_beadCount <- agrepl(noBB_sequences.df[row, 1], bead_sequence.df[seq, 2], max.distance = 1)
      antisense_beadCount <- agrepl(noBB_sequences.df[row, 2], bead_sequence.df[seq, 2], max.distance = 1)
      
      if (sense_matchCount){
        i = i + 1
      }
      if (antisense_matchCount){
        j = j + 1
      }
      if (sense_beadCount){
        p = p + 1
      }
      if (antisense_beadCount){
        q = q + 1
      }
    }
    # Append i, j, p, and q to noBB_sequences.df 
    noBB_sequences.df[row, 3] <- i
    noBB_sequences.df[row, 4] <- j
    noBB_sequences.df[row, 5] <- p
    noBB_sequences.df[row, 6] <- q
    
    # Add sense and antisense counts together
    for (row in 1:nrow(noBB_sequences.df)){
      noBB_sequences.df[row, 7] <- noBB_sequences.df[row, 3] + noBB_sequences.df[row, 4]
      noBB_sequences.df[row, 8] <- noBB_sequences.df[row, 5] + noBB_sequences.df[row, 6]
      
      # Change column names
      colnames(noBB_sequences.df)[3] <- "target_sense_count"
      colnames(noBB_sequences.df)[4] <- "target_antisense_count"
      colnames(noBB_sequences.df)[5] <- "bead_sense_count"
      colnames(noBB_sequences.df)[6] <- "bead_antisense_count"
      colnames(noBB_sequences.df)[7] <- "Target_Match_Sum"
      colnames(noBB_sequences.df)[8] <- "Bead_Match_Sum"
    }
    
    # End of sense and antisense partial match count function
  }
#---------------------------------------------------------------------------------------------------------------
# Print statements
  
cat(paste(crlf, crlf))
cat(paste("Sequences Grouped into Families With No Relation to Bead Sequences", crlf))

# Arrange noBB_sequences.df in descending order by copy number
noBB_sequences.df <- noBB_sequences.df[order(noBB_sequences.df[, 7], decreasing = TRUE), ]

# Loop to print each sequence on its own line
for (row in 1:nrow(noBB_sequences.df)){
  cat(paste("Sequence: ", noBB_sequences.df[row, 1], crlf, "Copies:", noBB_sequences.df[row, 7], crlf, "Matches in Bead Database: ", noBB_sequences.df[row, 8], crlf))
} 

# End of if dim(noBB.df > 0 if statement
}else{
  cat(paste("Sequences Grouped into Families With No Relation to Bead Sequences:", '\r'))
  cat("All protein sequences were found in families containing bead sequences.")
}
# End of if(dim(out_df[1] > 0 loop to make sure at least 1 family was found of the desired size))
}else{
  cat("No families of the desired size found.")
}


cat(crlf, crlf)
cat(paste("Alignment Complete"), crlf)

# Create dataframe for summary statistics
summary.stats.df <- data.frame("File Name" = character(0), "Total # of Sequences" = character(0), "# of N's" = character(0), "# of Exact Duplicates" = character(0), "# Number of close matches" = character(0), "Number of sense strands" = character(0), "Number of antisense strands" = character(0), "Number of shortmers" = character(0), "Number of longmers" = character(0), "Number Used in Alignment" = character(0), stringsAsFactors = FALSE)

# Load the dataframe with desired summary stats
summary.stats.df[1,1] <- alignment_fileName
summary.stats.df[1,2] <- allSeqs
summary.stats.df[1,3] <- number_of_N
summary.stats.df[1,4] <- no_of_dups
summary.stats.df[1,5] <- closeMatches
summary.stats.df[1,6] <- sense_before
summary.stats.df[1,7] <- anti_before
summary.stats.df[1,8] <- shortmers
summary.stats.df[1,9] <- longmers
summary.stats.df[1,10] <- no_for_alignment

cat(crlf, crlf)
cat(paste("Summary Statistics"), crlf)
cat(paste("Sequences in Original File  = "), allSeqs, crlf)
cat(paste("Sequences with at least 1 N = "), number_of_N, crlf)
cat(paste("Sense Strands = "), sense_before, crlf)
cat(paste("Antisense Strands = "), anti_before, crlf)
cat(paste("Exact Duplicates = "), no_of_dups, crlf)
cat(paste("Close Matches = "), closeMatches, crlf)
cat(paste("Shortmers (<37) = "), shortmers, crlf)
cat(paste("Longmers (>43) = "), longmers, crlf)
cat(paste("Aligned After Filtering = "), no_for_alignment, crlf)
# Print statement indicating that subsetting occurred automatically after all filters were applied to save time
if (subbed){
  cat(paste("Note: random subset of "), number_to_analyze, (" sequences used in the alignment."))
}
cat(crlf, crlf)

# Determine the most replicated sequences in the original sequence file (not from alignment, just from original file) and place in descending order
# Convert antisense to sense
replicate_frequency <- sequence.df
replicate_frequency.antisense <- sequence.df[agrep(fixed3, sequence.df$sequence, max.distance = 2, value = FALSE ),]
# Convert antisense sequences to their reverse complements
replicate_frequency.antisense$comp <- sapply(replicate_frequency.antisense$sequence, function(x) as.character(reverseComplement(DNAString(x))))
# Delete original sequences
replicate_frequency.antisense <- replicate_frequency.antisense[, -2]
# Rename antisense sequence column to match sense dataframe so antisense and sense dataframes can be merged
colnames(replicate_frequency.antisense)[2] <- "sequence"
# Subset sense sequences
replicate_frequency.sense <- sequence.df[agrep(fixed1, sequence.df$sequence, max.distance = 2, value = FALSE ),]
# Merge sense and antisense dataframes into a dataframe of sequences containing no perfect matches
replicate_frequency.merged <- bind_rows(replicate_frequency.sense, replicate_frequency.antisense)
replicates_df <- as.data.frame((table(replicate_frequency.merged$sequence)))
replicates_df <- replicates_df[order(replicates_df$Freq, decreasing = TRUE), ]
# Print top 8 most replicated sequences
cat(paste("Top 10 most replicated sequences:"), crlf)
for (i in 1:10){
  cat(paste(replicates_df[i, 1], ':', replicates_df[i, 2], crlf))
}


# Create csv file to export summary statistics (uncommenting this will overwrite the current csv file unless file name is changed!)
# write.csv(summary.stats.df,'/Users/danfeldheim/Documents/Data Files/iCX Alignments/AlignmentStats-proteins 9-104.csv', row.names = FALSE, col.names = TRUE)

# Append summary.stats.df to csv file
if(stat_append){
  write.table(summary.stats.df, file = stat_file, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  }

# Close the sink
sink()

# Create df for top ten sequence picks that emerged from alignment chosen by # of replicates
# Make sure they were non-bead binders by cross-referencing against original bead sequence file
if(dim(noBB.df)[1] > 0){
  
  # Create new df with first 10 rows (if noBB_sequences.df has 10 or more rows) and subset of columns from noBB.df
  # If noBB_sequences.df has fewer than 10 rows use all of them
  if(dim(noBB_sequences.df)[1] >= 10){
    sequence_picks.df <- noBB_sequences.df[1:10, c(1, 7, 8)]
  }
  else{
    sequence_picks.df <- noBB_sequences.df[, c(1, 7, 8)]
  }
  
  # Add column with file name
  sequence_picks.df$filename <- split_fileName
  # Reorder the columns
  sequence_picks.df <- sequence_picks.df[, c(4, 1, 2, 3)]
  # Add an empty row at the end of sequence_picks.df so a blank row prints between the top ten of each file
  sequence_picks.df[nrow(sequence_picks.df) + 1, ] <- ''

# Create csv file to export top 10 most replicated sequences that appear in families with nothing in common wtih bead database (uncommenting this will overwrite the current csv file!)
# write.csv(sequence_picks.df,'/Users/danfeldheim/Documents/Data Files/iCX Alignments/Sequence Picks.csv', row.names = FALSE, col.names = TRUE)

  if(append_picks){
# Append top 10 sequences to csv file
write.table(sequence_picks.df, file = pick_file, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  }
# End of if (dim(noBB.df)[1] > 0) from line 850
}

# Export text file of non bead binders (noBB.df[Row, 2])
# Open the sink to generate fasta text file from print statements below
# This will write a file that MEME can import, but R will not
if(nonbeadlist){
  nobeads <- paste0(report_path_and_file, '-nonbeadbinders')
  sink(nobeads, type = c("output", "message"))
  for (row in 1:nrow(noBB.df)){
    cat(paste0(">", noBB.df[row, 3], "\r"))
    cat(paste0(noBB.df[row, 2], "\r"))
  }
  # Close the sink
  sink()
 
}

end_time <- Sys.time()

run_time <- end_time - start_time
print (run_time)

# End of multiple-file processing for loop stared at top (line ~104)
}

closeAllConnections()




