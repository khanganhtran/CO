# setwd("~/Desktop/BINF/projects/CO/CO - SP")
# setwd("C:/Users/Ktran/Desktop/Bioinfo/CodonOptimization/Sanofi")
# install.packages("shiny", type ="binary") #This works

require(tidyr); require(dplyr); require(stringr); require(Biostrings); require(seqinr); require(openxlsx);require(ggplot2);require(e1071)

library(shiny)

# TO-DO
#[X] How to apply action button
#[X] rewrite code CodonOptimization() to take in partition (p) parameter
#[X] Output the CO table on RShiny
#[X] introduce a "save" button
#[ ] find out how to upload the app

#[X] Side Quest - download and link github

####################################################################################################


clean.seq <- function(seq){
  seq <- as.character(seq) #convert to character if not
  seq <- gsub("\\d*", "", seq) #remove all digits (if any)
  seq <- gsub("[[:space:]]", "", seq) #remove all white space, tabs, new line characters
  seq <- toupper(seq) #uppcase the aa string if not
  return(seq)
} #removes new line characters, numbers, and captilizes inputted string

read.fasta.DNA <- function(filename){
  filename <- as.character(filename)
  mySequences <- readDNAStringSet(filename) #Load the DataSet
  df <- data.frame(names(mySequences), paste(mySequences)); rm(mySequences)
  colnames(df) <- c("MRT", "sequence")
  df$MRT <- as.character(df$MRT); df$sequence <- as.character(df$sequence)
  
  return(df)
} #Reads in fasta seqs [DNA] and organizes into dataframe
read.fasta.AA <- function(filename){
  filename <- as.character(filename)
  mySequences <- readAAStringSet(filename) #Load the DataSet
  df <- data.frame(names(mySequences), paste(mySequences)); rm(mySequences)
  colnames(df) <- c("MRT", "sequence")
  df$MRT <- as.character(df$MRT); df$sequence <- as.character(df$sequence)
  
  return(df)
} #Reads in a fasta seqs [AA] and organizes into dataframe

split.AA <- function(seq,n) {
  
  #What if its 1 or 0 will through error
  
  if(n==0 | n==1){
    df <- data.frame(t(matrix(1:4)))
    colnames(df) <- c("n", "win.start", "win.end", "window")
    df$win.start <- 1
    df$win.end <- nchar(seq)
    df$window <- as.character(seq)
    return(df)
  }
  
  else{
    # Building df of coordinates
    coordinates <- 1:nchar(seq)
    L <- split(coordinates, cut(seq_along(coordinates), n, labels = FALSE))
    df <- data.frame(n = 1:length(L))
    
    df$win.start <- 0
    df$win.end <- 0
    df$window <- 0
    
    #use df of coordinates to extract sequence window(s)
    for (i in 1:length(L)) {
      df$win.start[i] <- L[[i]][1]
      
      x <- length(L[[i]])
      df$win.end[i] <- L[[i]][x]
      
      df$window[i] <- base::substring(seq, df$win.start[i], df$win.end[i])
    }
    
    return(df)
  }
}
DETECTOR <- function(sequence, motif, mismatch){
  seq <- sequence
  winsize <- nchar(motif)-1
  
  df <- seq(from = 1, to = nchar(seq), by = 1) #100
  df <- as.data.frame(df)
  df$x <- 0
  df$y <- 0
  df$seq <- 0
  df$mismatch <- 0
  for(i in 1:nchar(seq)){
    x <- i
    y <- i + winsize 
    frame <- substring(seq, x, y)
    
    if(nchar(frame) == winsize+1){
      df$x[i] <- x
      df$y[i] <- y
      df$seq[i] <- frame
    }
  }
  
  df<- df[df$seq != 0, c("x", "y", "seq", "mismatch")]
  motif.split <- strsplit(motif, "")[[1]]
  for (i in 1:nrow(df)) {
    frame.split <- strsplit(df$seq[i], "")[[1]]
    df$mismatch[i] <- hamming.distance(motif.split, frame.split) #Calculates Hamming distance between two seqs
  }
  
  df <- df[df$mismatch <= mismatch, ]
  
  if (dim(df)[1] == 0){
    return("No motif was located")
  }
  else{
    df <- df[order(df$mismatch), ]
    rownames(df) <- 1:nrow(df)
    
    df$pos <- 0
    
    for(i in 1:nrow(df)){
      df$pos[i] <- find.mm(seq1 = df$seq[i], seq2 = motif)
    }
    
    
    return(df)
  }
} #Detect desired Motif in sequence
find.mm <- function(seq1, seq2) {
  
  s1 <- as.character(seq1)
  s2 <- as.character(seq2)
  
  s1 <- as.data.frame(str_split(s1, pattern = ""))
  s2 <- as.data.frame(str_split(s2, pattern = ""))
  df.s <- data.frame(s1, s2)
  colnames(df.s) <- c("seq1", "seq2")
  df.s$pos <- 0
  
  for(i in 1:nrow(df.s)){
    if(df.s$seq1[i] != df.s$seq2[i]){
      df.s$pos[i] <- as.character(i)
    }
    else{
      df.s$pos[i] <- 0
    }
  }
  
  
  df.s <- df.s[df.s$pos != 0,]
  s.pos <- paste(df.s$pos, collapse = ",")
  
  return(s.pos)
  
}
GC.calc <- function(transcript) {
  sequence <- transcript
  G.count <- str_count(sequence, "G")
  C.count <- str_count(sequence, "C")
  seq.length <- nchar(sequence)
  
  answer <- (G.count + C.count) / seq.length
  return(answer)
} #GC calculator
GC.minmax <- function(string, size, minmax){
  
  #Bug in this part
  pat <- paste('(?<=.{', as.character(size), '})', sep = "")
  df <- as.data.frame(strsplit(string, pat, perl=TRUE))
  # return(df)
  df <- data.frame(1:nrow(df), df)
  colnames(df) <- c("n", "sequence")
  df$sequence <- as.character(df$sequence)
  
  df <- df %>% mutate(GC = GC.calc(sequence))
  
  
  
  minmax <- as.character(minmax)
  
  GC.vector <- as.vector(df$GC)
  GC.vector <- head(GC.vector, -1)
  
  if(minmax == "MAX"){
    # df <- df[1:nrow(df)-1,]
    return(max(GC.vector))
    # return("foo")
  }
  if(minmax == "MIN"){
    # df <- df[1:nrow(df)-1,]
    return(min(GC.vector))
  }
  else{
    return(GC.vector)
  }
}

splitInParts <- function(string, size){
  pat <- paste0('(?<=.{',size,'})')
  df <- as.data.frame(strsplit(string, pat, perl=TRUE))
  
  df <- data.frame(1:nrow(df), df)
  colnames(df) <- c("n", "sequence")
  df$sequence <- as.character(df$sequence)
  
  return(df)
}
CO <- function(aa, n, table){
  protein_seq <- as.character(aa)
  
  #Place AA string in a dataframe (df)
  n_trials <- as.integer(n)
  df <- tidyr::crossing(trial = 1:n_trials, seq_num = 1:nchar(protein_seq))
  df$aa = rep(unlist(str_split(protein_seq, "")), n_trials)
  
  df <- data.frame(as.data.frame(1:nrow(df)), df)
  colnames(df) <- c("n", "trial", "seq_num", "aa")
  
  #Step2: create pool of all amino amino acids (df.aa = for reading)
  df.human <- table
  # df.human <- read.xlsx("AminoAcidTable_Human_optimized_V2.xlsx", colNames = TRUE) #Load in AA chart
  df.aa <- as.data.frame(table(df$aa)) #table of amino acids to randomly generate
  colnames(df.aa) <- c("AA", "count")
  
  
  #start my initializing all dataframes
  df.stop <- 0; df.A <- 0; df.C <- 0 ; df.D <- 0; df.E <- 0; df.F <- 0; df.G <- 0; df.H <- 0; df.I <- 0; df.K <- 0; df.L <- 0; df.M <- 0; df.N <- 0; df.P <- 0; df.Q <- 0; df.R <- 0; df.S <- 0; df.T <- 0; df.V <- 0; df.W <- 0; df.Y <- 0
  
  ###Use 'Check' system to determine whether or not the make a dataframe
  
  #Check Whether not there are STOP codons
  check <- "*" %in%  unlist(str_split(protein_seq, ""))
  if (check == TRUE){
    df.stop <- sample(c(df.human$Codon[which(df.human$Codon == 'TAA')],
                        df.human$Codon[which(df.human$Codon == 'TAG')],
                        df.human$Codon[which(df.human$Codon == 'TGA')]), 
                      size = df.aa$count[df.aa$AA == "*"], prob = c(df.human$Frequency[which(df.human$Codon == 'TAA')],
                                                                    df.human$Frequency[which(df.human$Codon == 'TAG')],
                                                                    df.human$Frequency[which(df.human$Codon == 'TGA')]
                      ), replace = TRUE)
    df.stop <- cbind(df[df$aa == "*",], as.data.frame(df.stop))
    colnames(df.stop) <- c("n", "trial", "seq_num", "aa", "codon")
    df.stop$codon <- as.character(df.stop$codon)
  }
  check <- 0 #reset the counter
  
  #Check Whether not there are I (Isoleucine)
  check <- "I" %in%  unlist(str_split(protein_seq, ""))
  if (check == TRUE){ 
    df.I <- sample(c(df.human$Codon[which(df.human$Codon == 'ATT')],
                     df.human$Codon[which(df.human$Codon == 'ATC')],
                     df.human$Codon[which(df.human$Codon == 'ATA')]), 
                   size = df.aa$count[df.aa$AA == "I"], prob = c(df.human$Frequency[which(df.human$Codon == 'ATT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'ATC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'ATA')]
                   ), replace = TRUE)
    df.I <- cbind(df[df$aa == "I",], as.data.frame(df.I))
    colnames(df.I) <- c("n", "trial", "seq_num", "aa", "codon")
    df.I$codon <- as.character(df.I$codon)
  }
  check <- 0 #reset the counter
  
  #Check Whether not there are L (leucine)
  check <- "L" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.L <- sample(c(df.human$Codon[which(df.human$Codon == 'CTT')],
                     df.human$Codon[which(df.human$Codon == 'CTC')],
                     df.human$Codon[which(df.human$Codon == 'CTA')],
                     df.human$Codon[which(df.human$Codon == 'CTG')],
                     df.human$Codon[which(df.human$Codon == 'TTA')],
                     df.human$Codon[which(df.human$Codon == 'TTG')]), 
                   size = df.aa$count[df.aa$AA == "L"], prob = c(df.human$Frequency[which(df.human$Codon == 'CTT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CTC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CTA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CTG')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TTA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TTG')]
                   ), replace = TRUE)
    df.L <- cbind(df[df$aa == "L",], as.data.frame(df.L))
    colnames(df.L) <- c("n", "trial", "seq_num", "aa", "codon")
    df.L$codon <- as.character(df.L$codon)
  }
  check <- 0 #reset the counter
  
  
  #Check Whether not there are V (Valine)
  check <- "V" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.V <- sample(c(df.human$Codon[which(df.human$Codon == 'GTT')],
                     df.human$Codon[which(df.human$Codon == 'GTC')],
                     df.human$Codon[which(df.human$Codon == 'GTA')],
                     df.human$Codon[which(df.human$Codon == 'GTG')]), 
                   size = df.aa$count[df.aa$AA == "V"], prob = c(df.human$Frequency[which(df.human$Codon == 'GTT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GTC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GTA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GTG')]
                   ), replace = TRUE)
    df.V <- cbind(df[df$aa == "V",], as.data.frame(df.V))
    colnames(df.V) <- c("n", "trial", "seq_num", "aa", "codon")
    df.V$codon <- as.character(df.V$codon)
  }
  check <- 0 #reset the counter
  
  #Check Whether not there are F (Phenylalanine)
  check <- "F" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.F <- sample(c(df.human$Codon[which(df.human$Codon == 'TTT')],
                     df.human$Codon[which(df.human$Codon == 'TTC')]), 
                   size = df.aa$count[df.aa$AA == "F"], prob = c(df.human$Frequency[which(df.human$Codon == 'TTT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TTC')]
                   ), replace = TRUE)
    df.F <- cbind(df[df$aa == "F",], as.data.frame(df.F))
    colnames(df.F) <- c("n", "trial", "seq_num", "aa", "codon")
    df.F$codon <- as.character(df.F$codon)
  }
  check <- 0 #reset the counter
  
  
  #Check Whether not there are M (Methionine)
  check <- "M" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.M <- sample(c(df.human$Codon[which(df.human$Codon == 'ATG')]), 
                   size = df.aa$count[df.aa$AA == "M"], prob = c(df.human$Frequency[which(df.human$Codon == 'ATG')]
                   ), replace = TRUE)
    df.M <- cbind(df[df$aa == "M",], as.data.frame(df.M))
    colnames(df.M) <- c("n", "trial", "seq_num", "aa", "codon")
    df.M$codon <- as.character(df.M$codon)
  }
  check <- 0
  
  
  #Check Whether not there are C (Cysteine)
  check <- "C" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.C <- sample(c(df.human$Codon[which(df.human$Codon == 'TGT')],
                     df.human$Codon[which(df.human$Codon == 'TGC')]), 
                   size = df.aa$count[df.aa$AA == "C"], prob = c(df.human$Frequency[which(df.human$Codon == 'TGT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TGC')]
                   ), replace = TRUE)
    df.C <- cbind(df[df$aa == "C",], as.data.frame(df.C))
    colnames(df.C) <- c("n", "trial", "seq_num", "aa", "codon")
    df.C$codon <- as.character(df.C$codon)
  }
  check <- 0
  
  
  
  #Check Whether not there are A (Alanine)
  check <- "A" %in%  unlist(str_split(protein_seq, "")) 
  if(check == TRUE){
    df.A <- sample(c(df.human$Codon[which(df.human$Codon == 'GCT')],
                     df.human$Codon[which(df.human$Codon == 'GCC')],
                     df.human$Codon[which(df.human$Codon == 'GCA')],
                     df.human$Codon[which(df.human$Codon == 'GCG')]), 
                   size = df.aa$count[df.aa$AA == "A"], prob = c(df.human$Frequency[which(df.human$Codon == 'GCT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GCC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GCA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GCG')]
                   ), replace = TRUE)
    df.A <- cbind(df[df$aa == "A",], as.data.frame(df.A))
    colnames(df.A) <- c("n", "trial", "seq_num", "aa", "codon")
    df.A$codon <- as.character(df.A$codon)
  }
  check <- 0
  
  
  #Check Whether not there are G (Glycine)
  check <- "G" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.G <- sample(c(df.human$Codon[which(df.human$Codon == 'GGT')],
                     df.human$Codon[which(df.human$Codon == 'GGC')],
                     df.human$Codon[which(df.human$Codon == 'GGA')],
                     df.human$Codon[which(df.human$Codon == 'GGG')]), 
                   size = df.aa$count[df.aa$AA == "G"], prob = c(df.human$Frequency[which(df.human$Codon == 'GGT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GGC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GGA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GGG')]
                   ), replace = TRUE)
    df.G <- cbind(df[df$aa == "G",], as.data.frame(df.G))
    colnames(df.G) <- c("n", "trial", "seq_num", "aa", "codon")
    df.G$codon <- as.character(df.G$codon)
  }
  check <- 0
  
  #Check Whether not there are P (Proline)
  check <- "P" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.P <- sample(c(df.human$Codon[which(df.human$Codon == 'CCT')],
                     df.human$Codon[which(df.human$Codon == 'CCC')],
                     df.human$Codon[which(df.human$Codon == 'CCA')],
                     df.human$Codon[which(df.human$Codon == 'CCG')]), 
                   size = df.aa$count[df.aa$AA == "P"], prob = c(df.human$Frequency[which(df.human$Codon == 'CCT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CCC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CCA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CCG')]
                   ), replace = TRUE)
    df.P <- cbind(df[df$aa == "P",], as.data.frame(df.P))
    colnames(df.P) <- c("n", "trial", "seq_num", "aa", "codon")
    df.P$codon <- as.character(df.P$codon)
  }
  check <- 0
  
  #Check Whether not there are T (Threonine)
  check <- "T" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.T <- sample(c(df.human$Codon[which(df.human$Codon == 'ACT')],
                     df.human$Codon[which(df.human$Codon == 'ACC')],
                     df.human$Codon[which(df.human$Codon == 'ACA')],
                     df.human$Codon[which(df.human$Codon == 'ACG')]), 
                   size = df.aa$count[df.aa$AA == "T"], prob = c(df.human$Frequency[which(df.human$Codon == 'ACT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'ACC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'ACA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'ACG')]
                   ), replace = TRUE)
    df.T <- cbind(df[df$aa == "T",], as.data.frame(df.T))
    colnames(df.T) <- c("n", "trial", "seq_num", "aa", "codon")
    df.T$codon <- as.character(df.T$codon)
  }
  check <- 0
  
  #Check Whether not there are S (Serine)
  check <- "S" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.S <- sample(c(df.human$Codon[which(df.human$Codon == 'TCT')],
                     df.human$Codon[which(df.human$Codon == 'TCC')],
                     df.human$Codon[which(df.human$Codon == 'TCA')],
                     df.human$Codon[which(df.human$Codon == 'TCG')],
                     df.human$Codon[which(df.human$Codon == 'AGT')],
                     df.human$Codon[which(df.human$Codon == 'AGC')]), 
                   size = df.aa$count[df.aa$AA == "S"], prob = c(df.human$Frequency[which(df.human$Codon == 'TCT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TCC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TCA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TCG')],
                                                                 df.human$Frequency[which(df.human$Codon == 'AGT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'AGC')]
                   ), replace = TRUE)
    df.S <- cbind(df[df$aa == "S",], as.data.frame(df.S))
    colnames(df.S) <- c("n", "trial", "seq_num", "aa", "codon")
    df.S$codon <- as.character(df.S$codon)
  }
  check <- 0
  
  
  #Check Whether not there are Y (Tyrosine)
  check <- "Y" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.Y <- sample(c(df.human$Codon[which(df.human$Codon == 'TAT')],
                     df.human$Codon[which(df.human$Codon == 'TAC')]), 
                   size = df.aa$count[df.aa$AA == "Y"], prob = c(df.human$Frequency[which(df.human$Codon == 'TAT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'TAC')]
                   ), replace = TRUE)
    df.Y <- cbind(df[df$aa == "Y",], as.data.frame(df.Y))
    colnames(df.Y) <- c("n", "trial", "seq_num", "aa", "codon")
    df.Y$codon <- as.character(df.Y$codon)
  }
  check <- 0
  
  
  #Check Whether not there are W (Tryptophan)
  check <- "W" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.W <- sample(c(df.human$Codon[which(df.human$Codon == 'TGG')]), 
                   size = df.aa$count[df.aa$AA == "W"], prob = c(df.human$Frequency[which(df.human$Codon == 'TGG')]
                   ), replace = TRUE)
    df.W <- cbind(df[df$aa == "W",], as.data.frame(df.W))
    colnames(df.W) <- c("n", "trial", "seq_num", "aa", "codon")
    df.W$codon <- as.character(df.W$codon)
  }
  check <- 0
  
  
  #Check Whether not there are Q (Glutamine)
  check <- "Q" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.Q <- sample(c(df.human$Codon[which(df.human$Codon == 'CAA')],
                     df.human$Codon[which(df.human$Codon == 'CAG')]), 
                   size = df.aa$count[df.aa$AA == "Q"], prob = c(df.human$Frequency[which(df.human$Codon == 'CAA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CAG')]
                   ), replace = TRUE)
    df.Q <- cbind(df[df$aa == "Q",], as.data.frame(df.Q))
    colnames(df.Q) <- c("n", "trial", "seq_num", "aa", "codon")
    df.Q$codon <- as.character(df.Q$codon)
  }
  check <- 0
  
  
  #Check Whether not there are N (Asparagine)
  check <- "N" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.N <- sample(c(df.human$Codon[which(df.human$Codon == 'AAT')],
                     df.human$Codon[which(df.human$Codon == 'AAC')]), 
                   size = df.aa$count[df.aa$AA == "N"], prob = c(df.human$Frequency[which(df.human$Codon == 'AAT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'AAC')]
                   ), replace = TRUE)
    df.N <- cbind(df[df$aa == "N",], as.data.frame(df.N))
    colnames(df.N) <- c("n", "trial", "seq_num", "aa", "codon")
    df.N$codon <- as.character(df.N$codon)
  }
  check <- 0
  
  #Check Whether not there are H (Histidine)
  check <- "H" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.H <- sample(c(df.human$Codon[which(df.human$Codon == 'CAT')],
                     df.human$Codon[which(df.human$Codon == 'CAC')]), 
                   size = df.aa$count[df.aa$AA == "H"], prob = c(df.human$Frequency[which(df.human$Codon == 'CAT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CAC')]
                   ), replace = TRUE)
    df.H <- cbind(df[df$aa == "H",], as.data.frame(df.H))
    colnames(df.H) <- c("n", "trial", "seq_num", "aa", "codon")
    df.H$codon <- as.character(df.H$codon)
  }
  check <- 0
  
  #Check Whether not there are E (Glutamic Acid)
  check <- "E" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.E <- sample(c(df.human$Codon[which(df.human$Codon == 'GAA')],
                     df.human$Codon[which(df.human$Codon == 'GAG')]), 
                   size = df.aa$count[df.aa$AA == "E"], prob = c(df.human$Frequency[which(df.human$Codon == 'GAA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GAG')]
                   ), replace = TRUE)
    df.E <- cbind(df[df$aa == "E",], as.data.frame(df.E))
    colnames(df.E) <- c("n", "trial", "seq_num", "aa", "codon")
    df.E$codon <- as.character(df.E$codon)
  }
  check <- 0
  
  
  #Check Whether not there are D (Aspartic Acid)
  check <- "D" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.D <- sample(c(df.human$Codon[which(df.human$Codon == 'GAT')],
                     df.human$Codon[which(df.human$Codon == 'GAC')]), 
                   size = df.aa$count[df.aa$AA == "D"], prob = c(df.human$Frequency[which(df.human$Codon == 'GAT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'GAC')]
                   ), replace = TRUE)
    df.D <- cbind(df[df$aa == "D",], as.data.frame(df.D))
    colnames(df.D) <- c("n", "trial", "seq_num", "aa", "codon")
    df.D$codon <- as.character(df.D$codon)
  }
  check <- 0
  
  #Check Whether not there are K (Lysine)
  check <- "K" %in%  unlist(str_split(protein_seq, ""))
  if(check == TRUE){
    df.K <- sample(c(df.human$Codon[which(df.human$Codon == 'AAA')],
                     df.human$Codon[which(df.human$Codon == 'AAG')]), 
                   size = df.aa$count[df.aa$AA == "K"], prob = c(df.human$Frequency[which(df.human$Codon == 'AAA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'AAG')]
                   ), replace = TRUE)
    df.K <- cbind(df[df$aa == "K",], as.data.frame(df.K))
    colnames(df.K) <- c("n", "trial", "seq_num", "aa", "codon")
    df.K$codon <- as.character(df.K$codon)
  }
  check <- 0
  
  
  
  #Check Whether not there are R (Arginine)
  check <- "R" %in%  unlist(str_split(protein_seq, ""))
  
  if(check == TRUE){
    df.R <- sample(c(df.human$Codon[which(df.human$Codon == 'CGT')],
                     df.human$Codon[which(df.human$Codon == 'CGC')],
                     df.human$Codon[which(df.human$Codon == 'CGA')],
                     df.human$Codon[which(df.human$Codon == 'CGG')],
                     df.human$Codon[which(df.human$Codon == 'AGA')],
                     df.human$Codon[which(df.human$Codon == 'AGG')]), 
                   size = df.aa$count[df.aa$AA == "R"], prob = c(df.human$Frequency[which(df.human$Codon == 'CGT')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CGC')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CGA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'CGG')],
                                                                 df.human$Frequency[which(df.human$Codon == 'AGA')],
                                                                 df.human$Frequency[which(df.human$Codon == 'AGG')]
                   ), replace = TRUE)
    df.R <- cbind(df[df$aa == "R",], as.data.frame(df.R))
    colnames(df.R) <- c("n", "trial", "seq_num", "aa", "codon")
    df.R$codon <- as.character(df.R$codon)
  }
  check <- 0
  
  
  df <- rbind(df.stop, df.A, df.C, df.D, df.E, df.F, df.G, df.H, df.I, df.K, df.L, df.M, df.N, df.P, df.Q, df.R, df.S, df.T, df.V, df.W, df.Y)
  df <- df[order(df$n), ]
  df <- df[df$n != 0,]
  
  
  #How to group a dataframe by a trial and then collapse the Codon into a string
  df <- (df %>% 
           group_by(trial) %>%
           summarize(mrna = paste(codon, collapse="")))
  
  
  rm(df.aa, df.stop, df.A, df.C, df.D, df.E, df.F, df.G, df.H, df.I, df.K, df.L, df.M, df.N, df.P, df.Q, df.R, df.S, df.T, df.V, df.W, df.Y)
  
  return(df)
} #inputs aa string, n desired number of seq, freq table. Outputs df (n, mRNA)
CO.filter <- function(seq, n, CAI, codons, nCIS){
  
  seq <- as.character(seq)
  CAI <- as.numeric(CAI)
  n <- as.integer(n)
  
  df.human <- as.character(codons)
  df.nCIS <- as.character(nCIS)
  
  df.human <- read.xlsx(df.human, colNames = TRUE) #Load in AA chart
  w <- df.human[order(df.human$Codon),] #Preparing vector of RelativeAdaptiveness to calculate CAI
  rownames(w) <- tolower(w$Codon)
  w <- w$RelativeAdaptiveness 
  
  df.nCIS <- read.fasta.DNA(df.nCIS) #Load in nCIS elements
  
  
  df <- CO(aa = seq, n = n, table = df.human) #Codon optimize n sequences
  
  df <- df %>% mutate(GC = GC.calc(mrna)) #calculate transcript GC
  
  #Establishing analytical variables
  df[,c("CAI", "CFD", "GC.max", "GC.min", "motif.ATCTGTT", "motif.TTTTTT", "motif.HINDI", "motif.Sap1", "motif.Sap12", "motif.Xbal", "motif.aarl1", "motif.aarl2")] <- 0
  
  #Calculating transcript analytics
  for(i in 1:nrow(df)){
    
    #Calculating GC max/min in windows
    GC.max <- GC.minmax(string = df$mrna[i], size = 30, minmax = "MAX")
    df$GC.max[i] <- GC.max
    
    GC.min <- GC.minmax(string = df$mrna[i], size = 30, minmax = "MIN")
    df$GC.min[i] <- GC.min
    
    ####### In-house Motif analysis
    object <- DNAString(df$mrna[i]) #DNAString object to be used for motif analysis
    
    #rrnB1 terminator
    df$motif.ATCTGTT[i] <- m.motif(object, motif = "ATCTGTT", max = 1, min = 0)
    
    #Stretch of T's
    df$motif.TTTTTT[i] <- m.motif(object, motif = "TTTTTT", max = 0, min = 0)
    
    #HindIII - AAGCTT
    df$motif.HINDI[i] <- m.motif(object, motif = "AAGCTT", max = 0, min = 0)
    
    #Sap1 - GAAGAGC
    df$motif.Sap1[i] <- m.motif(object, motif = "GAAGAGC", max = 0, min = 0)
    df$motif.Sap12[i] <- m.motif(object, motif = "GCTCTTC", max = 0, min = 0)
    
    #Xbal - TCTAGA
    df$motif.Xbal[i] <- m.motif(object, motif = "TCTAGA", max = 0, min = 0)
    
    #AarI (CACCTGC and GCAGGTG, exact match)
    df$motif.aarl1[i] <- m.motif(object, motif = "CACCTGC", max = 0, min = 0)
    df$motif.aarl2[i] <- m.motif(object, motif = "GCAGGTG", max = 0, min = 0)
    
    #Analyzing CAI
    mrna.parse <- str_split(tolower(df$mrna[i]), "")[[1]]
    df$CAI[i] <- cai(mrna.parse, w, numcode = 1, zero.threshold = 0.0001, zero.to = 0.01)
    
    #Analyzing CFD
    df$CFD[i] <- assign.CFD(seq = df$mrna[i], table = df.human)
    
    rm(object) #remove object for aesthetic purposes #put this OUTSIDE
    
  }
  
  #Filtering the transcripts (analytics)
  df <- df[df$GC.max < 0.70 &
             df$GC.min > 0.30 &
             df$motif.ATCTGTT == 0 &
             df$motif.TTTTTT == 0 &
             df$motif.HINDI == 0 &
             df$motif.Sap1 == 0 &
             df$motif.Sap12 == 0 &
             df$motif.Xbal == 0 &
             
             df$motif.aarl1 == 0 &
             df$motif.aarl2 == 0 &
             
             df$CAI >= CAI,] #Filter out the "bad" transcripts
  
  df<- unique(df) #grab only UNIQUE transcripts!
  df <- df[order(-df$CAI),]
  
  df <- find.nCIS(df, df.nCIS) #analyze inputed df for negative CIS elements
  
  df$nRE <- 0 #analyze df for negative Repeat elements
  for(i in 1:nrow(df)){
    object <- DNAString(df$mrna[i]) #DNAString object to be used for motif analysis
    
    df$nRE[i] <- count.NRE(sequence = df$mrna[i])
    
    rm(object)
  }
  
  return(df)
  
} #inputs aa string, n desired number of seq, freq table. Outputs df (n, mRNA, GC, CAI, CFD, nCIS, nRE)
codonoptimize <- function(seq, p, n, CAI, codons, nCIS){
  
  start.time <- Sys.time() #Start the time
  
  seq <- clean.seq(seq)
  p <- as.integer(p)
  CAI <- as.numeric(CAI)
  n <- as.integer(n)
  
  df.human <- as.character(codons)
  df.nCIS <- as.character(nCIS)
  
  print("Codon Optimizing Windows...")
  
  df.prot <- split.AA(seq=seq, n = p) 
  df.master <- NULL
  
  initial <- 0.2
  
  withProgress(message = 'Codon Optimizing Windows...', style = "notification", detail = "window 1", value = initial, {
    for (i in 1:nrow(df.prot)) {
      
      df <- CO.filter(seq = df.prot$window[i], 
                      n=n, 
                      CAI=CAI, 
                      codons=df.human, 
                      nCIS=df.nCIS
      )
      #if df.master is empty, it becomes the first df
      if(isEmpty(df.master[1]) == TRUE){
        df.master <- df[,c("n", "mrna")]
        #Write code there that will clean up the df
        colnames(df.master) <- c("n", "mrna")
      }
      else{
        # if not, append subsequent dataframes
        df <- df[,c("n", "mrna")]
        df.master <- merge(df.master, df, by = "n", all.y = TRUE, all.x = TRUE)
        colnames(df.master) <- c("n", sprintf("window%02d", seq(1,ncol(df.master)-1)))
      }
      
      text <- paste("Window", as.character(i), "finished!", sep = " ")
      print(text)
      withProgress(message = text,
                   style = "notification", value = NULL, {
                     Sys.sleep(1) #This suspends the execution of futher code for n seconds
                   })
  

      # incProgress(1/nrow(df.prot), detail = paste("window", i+1))
      
      incProgress((1/nrow(df.prot))*(1-initial), detail = paste("Window", i+1))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.25)
    }
  })
  
  ############################################################
  
  print("All windows have finished!")
  
  time <- as.character(round(Sys.time()-start.time, 2))
  print(paste("Time elapsed:", time, sep = " " ))
  
  
  withProgress(message = 'All windows have finished!',
               style = "notification", detail = paste("Time elapsed:", time, sep = " " ), value = NULL, {
                 Sys.sleep(2) #This suspends the execution of futher code for n seconds
               })
  
  
  

  
  df.human <- read.xlsx(df.human, colNames = TRUE) #Load in AA chart
  w <- df.human[order(df.human$Codon),] #Preparing vector of RelativeAdaptiveness to calculate CAI
  rownames(w) <- tolower(w$Codon)
  w <- w$RelativeAdaptiveness 
  
  df.nCIS <- read.fasta.DNA(df.nCIS) #Load in nCIS elements
  
  rm(df)
  rm(df.prot)
  
  
  df.master <- df.master[complete.cases(df.master), ] #Remove dataframes with na values
  
  
  #Create random dataframe with n rows
  df.apprentice <- as.data.frame(matrix(0, ncol = ncol(df.master)-1, nrow = 300))
  #Populate dataframe with integers in df.master (biased towards higher half)
  for(i in 1:ncol(df.apprentice)){
    df.apprentice[,i] <- sample(x = floor(nrow(df.master)/2), size = 50, replace = TRUE)
  }
  
  #Double for-loop to use coordinates from df.app to grab seq-windows from df.master
  for(i in 1:nrow(df.apprentice)){
    for(j in 1:ncol(df.apprentice)){
      n <- df.apprentice[i,j] #the number at the incremental coordinates
      df.apprentice[i,j] <- df.master[,-1][n, j]
    }
  }
  
  #Combine df.apprentice & df.master
  colnames(df.apprentice) <- c(sprintf("window%02d", seq(1,ncol(df.master)-1)))
  df.master <- rbind(df.master[,-1] , df.apprentice)
  
  rm(df.apprentice)
  
  df.master$mrna <- 0
  for(i in 1:nrow(df.master)){
    seq <- paste(df.master[i,1:ncol(df.master)-1], collapse = "")
    df.master$mrna[i] <- seq
  }
  
  #Final df of CO sequences (Final analytics need to be performed)
  df.master <- data.frame(n = 1:nrow(df.master), mrna = df.master$mrna)
  df.master <- data.frame(n = 1:length(unique(df.master$mrna)), mrna = unique(df.master$mrna)) #only unique candidates
  
  df <- df.master; rm(df.master)
  df$mrna <- as.character(df$mrna)
  
  print("Table of Codon Optimized mRNA made!")
  
  withProgress(message = 'Table of Codon Optimized mRNA made!',
               style = "notification", detail = "Preparing Analytics", value = NULL, {
                 Sys.sleep(2) #This suspends the execution of futher code for n seconds
               })
  

  
  
  df <- df %>% mutate(GC = GC.calc(mrna)) #calculate transcript GC
  
  #Establishing analytical variables
  df[,c("CAI", "CFD", "GC.max", "GC.min", "motif.ATCTGTT", "motif.TTTTTT", "motif.HINDI", "motif.Sap1", "motif.Sap12", "motif.Xbal", "motif.aarl1", "motif.aarl2")] <- 0
  
  print("Performing mRNA analytics...")
  
  withProgress(message = 'Calculating mRNA analytics...',
               style = "notification", detail = "Preparing Analytics", value = NULL, {
                 Sys.sleep(2) #This suspends the execution of futher code for n seconds
  #Calculating transcript analytics
  for(i in 1:nrow(df)){
    #Calculating GC max/min in windows
    GC.max <- GC.minmax(string = df$mrna[i], size = 30, minmax = "MAX")
    df$GC.max[i] <- GC.max
    GC.min <- GC.minmax(string = df$mrna[i], size = 30, minmax = "MIN")
    df$GC.min[i] <- GC.min
    object <- DNAString(df$mrna[i]) #DNAString object to be used for motif analysis
    #rrnB1 terminator
    df$motif.ATCTGTT[i] <- m.motif(object, motif = "ATCTGTT", max = 1, min = 0)
    #Stretch of T's
    df$motif.TTTTTT[i] <- m.motif(object, motif = "TTTTTT", max = 0, min = 0)
    #HindIII - AAGCTT
    df$motif.HINDI[i] <- m.motif(object, motif = "AAGCTT", max = 0, min = 0)
    #Sap1 - GAAGAGC
    df$motif.Sap1[i] <- m.motif(object, motif = "GAAGAGC", max = 0, min = 0)
    df$motif.Sap12[i] <- m.motif(object, motif = "GCTCTTC", max = 0, min = 0)
    #Xbal - TCTAGA
    df$motif.Xbal[i] <- m.motif(object, motif = "TCTAGA", max = 0, min = 0)
    #AarI (CACCTGC and GCAGGTG, exact match)
    df$motif.aarl1[i] <- m.motif(object, motif = "CACCTGC", max = 0, min = 0)
    df$motif.aarl2[i] <- m.motif(object, motif = "GCAGGTG", max = 0, min = 0)
    #Analyzing CAI
    mrna.parse <- str_split(tolower(df$mrna[i]), "")[[1]]
    df$CAI[i] <- cai(mrna.parse , w, numcode = 1, zero.threshold = 0.0001, zero.to = 0.01)
    #Analyzing CFD
    df$CFD[i] <- assign.CFD(seq = df$mrna[i], table = df.human)
    
    rm(object) #remove object for aesthetic purposes #put this OUTSIDE
    percent <- as.character(round(i/nrow(df)*100),5)
    incProgress(1/nrow(df), detail = paste(percent, "%")) #Tricky
    Sys.sleep(1/nrow(df))
  }
               })
  
  withProgress(message = 'Calculating mRNA analytics...',
               style = "notification", detail = "Filtering Transcripts", value = NULL, {
                 Sys.sleep(2) #This suspends the execution of futher code for n seconds
               })
  
  #Filtering the transcripts (analytics)
  df <- df[
    df$GC.max < 0.70 &
      df$GC.min > 0.30 &
      df$motif.ATCTGTT == 0 &
      df$motif.TTTTTT == 0 &
      df$motif.HINDI == 0 &
      df$motif.Sap1 == 0 &
      df$motif.Sap12 == 0 &
      df$motif.Xbal == 0 &
      
      df$motif.aarl1 == 0 &
      df$motif.aarl2 == 0 &
      
      df$CAI >= CAI,] #Filter out the "bad" transcripts
  
  df<- unique(df) #grab only UNIQUE transcripts!
  df <- df[order(-df$CAI),]
  
  df <- find.nCIS(df, df.nCIS) #analyze inputed df for negative CIS elements
  
  withProgress(message = 'Calculating mRNA analytics...',
               style = "notification", detail = "Counting Negative Regulatory elements", value = NULL, {
                 Sys.sleep(2) #This suspends the execution of futher code for n seconds
               })
  

  
  df$nRE <- 0 #analyze df for negative Repeat elements
  for(i in 1:nrow(df)){
    object <- DNAString(df$mrna[i]) #DNAString object to be used for motif analysis
    
    df$nRE[i] <- count.NRE(sequence = df$mrna[i])
    
    rm(object)
  }
  
  withProgress(message = 'Calculating mRNA analytics...',
               style = "notification", detail = "Counting Negative Repeats", value = NULL, {
                 Sys.sleep(2) #This suspends the execution of futher code for n seconds
               })
  
  
  print("Transcript Analytics Finished!")
  
  time <- as.character(round(Sys.time()-start.time, 2))
  print(paste("Total Run-time:", time, sep = " " ))
  
  
  withProgress(message = 'Transcript Analytics Finished!',
               style = "notification", detail = paste("Total Run-time:", time, sep = " " ), value = NULL, {
                 Sys.sleep(4) #This suspends the execution of futher code for n seconds
               })
  
  
  rm(list = ls.str(mode = 'numeric')) #remove all numeric variables
  rm(list = ls.str(mode = 'character')) #remove all character variables
  
  return(df)
  
  
}

assign.rarity <- function(codon, df){
  
  codon <- as.character(codon)
  df.rare <- df
  
  #Check parameter
  if(codon %in% df.rare$Codon == TRUE){
    return(df.rare[df.rare$Codon == codon, ][["Status"]])
  }
  
} #Function assigns "rarity" of inputted codon based on freq table
assign.CFD <- function(seq, table){
  seq <- as.character(seq)
  
  #Parse the seq into codon variables
  seq.df <- as.data.frame(strsplit(seq, '(?<=.{3})', perl=TRUE))
  colnames(seq.df) <- "Codon"
  seq.df$Codon <- as.character(seq.df$Codon)
  
  #initialize rarity variable
  seq.df$Status <- 0
  for(i in 1:nrow(seq.df)){
    #Assign rarity status to each variable
    seq.df$Status[i] <- assign.rarity(seq.df$Codon[i], table)
  }
  
  seq.df$Status <- gsub("rare", "R", seq.df$Status)
  seq.df$Status <- gsub("common", "c", seq.df$Status)
  
  #seq to be analyzed for CFD
  seq.rare <- paste(seq.df$Status, collapse = "")
  
  #calculate CFD (#Rare codon duplets / AA string)
  CFD.RR <- nrow(as.data.frame(matchPattern("RR", DNAString(seq.rare), max.mismatch = 0, min.mismatch = 0)))
  CFD.RcR <- nrow(as.data.frame(matchPattern("RcR", DNAString(seq.rare), max.mismatch = 0, min.mismatch = 0)))
  
  CFD <- (CFD.RR + CFD.RcR) / (nchar(seq)/3)
  
  return(round(CFD, 3))
  
} #Calculates CFD from inputed mRNA based on freq table
m.motif <- function(seq, motif, max, min){
  seq <- seq
  motif <- as.character(motif)
  
  if(mode(seq) != "S4"){
    #If it is not a DNAString obj --> make it one
    seq <- Biostrings::DNAString(seq)
  }
  
  #Find motif and return #
  df <- as.data.frame(matchPattern(motif, seq, max.mismatch = max, min.mismatch = min))
  return(nrow(df))
} #Function for motif (object/string, motif, max, min)
find.nCIS <- function(df, df2){
  
  df.nCIS <- df2
  
  #Grab only the CO sequences & relevant analytics (GC, CAI, CFD)
  df.temp <- data.frame(1:nrow(df), df[,c("mrna", "GC", "CAI", "CFD")]) #inputted df must have mrna, GC, CAI, CFD
  colnames(df.temp)[1] <- "n"
  
  
  
  #I think we're going to have code by analyzing each CIS motif against filtered transcripts
  for(i in 1:nrow(df.nCIS)){
    motif <- df.nCIS$sequence[i]
    motif <- as.data.frame(grepl(motif, df$mrna))
    df.temp <- cbind(df.temp, motif)
  } 
  
  
  df.temp2 <- df.temp[,6:ncol(df.temp)] #df of only TRUE/FALSE CIS motifs
  df.temp <- df.temp[,1:5]#df of only mRNA + other analytics
  
  # colnames(df.temp2) <- sprintf("M%02d", seq(1,ncol(df.temp2))) #we might not need this line
  
  df.temp2 <- data.frame(apply(df.temp2, 1, function(x) length(which(x== TRUE)))) #This works
  colnames(df.temp2) <- "nCIS"
  
  df.temp <- data.frame(df.temp, df.temp2) #a success!!!!!
  
  return(df.temp)
  
} #analyze inputed df for negative CIS elements

#Functions to find NegRepeatElements
findTandemRepeats0 <- function(subject, period.multiple=2, include.period1=FALSE, min.length=24) {
  if (!isSingleNumber(period.multiple))
    stop("'period.multiple' must be a single integer")
  if (!is.integer(period.multiple))
    period.multiple <- as.integer(period.multiple)
  if (period.multiple < 2L)
    stop("'period.multiple' must be >= 2")
  if (!isTRUEorFALSE(include.period1))
    stop("'include.period1' must be TRUE or FALSE")
  if (!isSingleNumber(min.length))
    stop("'min.length' must be a single integer")
  if (!is.integer(min.length))
    min.length <- as.integer(min.length)
  if (min.length < 12L)
    stop("'min.length' must be >= 12")
  s1 <- subseq(subject, start=1L+period.multiple)
  s2 <- subseq(subject, end=-1L-period.multiple)
  ir <- as(as.raw(s1) == as.raw(s2), "NormalIRanges")
  ir <- ir[width(ir) >= period.multiple]
  end(ir) <- end(ir) + period.multiple
  ir <- ir[width(ir) >= min.length]
  repeats <- Views(subject, ir)
  
  af <- alphabetFrequency(repeats, baseOnly=TRUE)
  ok <- af[ , "other"] == 0L  # has no IUPAC ambiguities
  if (!include.period1)
    ok <- ok & rowSums(af[ , DNA_BASES] != 0L) >= 2L
  repeats[which(ok)]
}
findTandemRepeats <- function(subject, include.period1=FALSE) {
  trs_list <- lapply(7:12,
                     function(period.multiple)
                       ranges(findTandemRepeats0(subject,
                                                 period.multiple,
                                                 include.period1))
  )
  trs <- sort(unique(do.call("c", trs_list)))
  Views(subject, trs)
} ## Find all *exact* tandem repeats of period <= 24.
count.NRE <- function(sequence){
  #Limit: minimum input sequence 10, NRE must be 24 nt minimum 
  sequence <- sequence
  
  if(mode(sequence) != "S4"){
    #If it is not a DNAString obj --> make it one
    sequence <- Biostrings::DNAString(sequence)
  }
  
  n <- lengths(as.vector(as.data.frame(findTandemRepeats(sequence, include.period1=TRUE))))
  return(n)
  
} #input sequence + return NRE count (if any); 
find.NRE <- function(sequence){
  sequence <- sequence
  
  if(mode(sequence) != "S4"){
    #If it is not a DNAString obj --> make it one
    sequence <- Biostrings::DNAString(sequence)
  }
  
  n <- findTandemRepeats(sequence, include.period1=TRUE)
  return(n)
} #Input sequence + return NRE (if any)


translate.seq <- function(seq){
  seq <- as.character(Biostrings::translate(DNAString(as.character(seq))))
  return(seq)
} #translate DNA sequence


####################################################################################################

# TO-DO
#[X] How to apply action button
#rewrite code CodonOptimization() to take in partition (p) parameter


ui <- fluidPage(
  # titlePanel("TBio Codon Optimization"),

  titlePanel(title =  div(img(src="logo.jpg", height="25%", width="25%"))),
  
  # titlePanel(fluidRow(
  #   column(9, "Codon Optimization"), 
  #   column(3, img(height="50%", width="50%", src = "logo.jpg"))
  # )
  # ),
  
  
  sidebarLayout(
    sidebarPanel(
      h4(strong("Codon Optimization tool"), align = "center"),
      
      # helpText("Input codon optimization parameters","(sequence should be captilized)"),
      
      textInput("aa", label = "Input amino acid sequence", value = "MATGSRTSLLLAFGLLCLPWLQEGSAFPTIPLSQSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQSEDEADYYCNSLTSISTWVFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS*" ),
      selectInput("p", label = "Sequence partition(s)", choices = as.character(c(1:10)), selected = "#"),
      selectInput("n", label = "Generate sequence(s)", choices = as.character(seq(1000,10000, by = 1000)), selected = "#"),
      selectInput("cai", label = "Codon Adaptation Index (min)", choices = as.character(seq(0.6,0.85, by = 0.05)), selected = "#"),
      
      # selectInput("cai", label = "Sequence Partition(s)", choices = as.character(c(1:10)), selected = "#"),
      
      fileInput("codon", label = "Codon Frequency Table", placeholder = "*.xlsx", accept = c(".xlsx")),
      fileInput("nCIS", label = "Negative CIS file", placeholder = "*.fasta"),
      
      br(),
      actionButton("go", "submit"),
      
      downloadButton("downloadData", "Download"),
      # actionButton("go2", "submit"),
      # radioButtons('style', 'Progress bar style', c('notification')),
      helpText("Click 'submit' to codon optimize amino acid sequence")
      ),
    
    mainPanel(
      # textOutput("AA"), #amino acid String
      # textOutput("P"), #Partition Number
      # 
      # textOutput("N"), #generated sequences
      # 
      # textOutput("CAI"),
      # 
      # tableOutput("CODON"),
      # tableOutput("NCIS"),
      
      tableOutput("CO") #the holy Grail

    )
  )
  )

server <- function(input, output) {
  #Save the variables (as eventReactive)
  aa <- eventReactive(input$go, {input$aa})
  p <- eventReactive(input$go, {input$p})
  n <- eventReactive(input$go, {input$n})
  cai <- eventReactive(input$go, {input$cai})
  codon <- eventReactive(input$go, {input$codon})
  NCIS <- eventReactive(input$go, {input$nCIS})
  
  #Specifying arguments
  # output$AA <- renderText({paste("You have selected:", aa())})
  # output$P <- renderText({paste("Number of Partions:", p())})
  # output$N <- renderText({paste("Number of sequences:", n())})
  # output$CAI <- renderText({paste("Minimum CAI:", cai())})
  # 
  # output$CODON <- renderTable({
  #   req(codon())
  #   inFile <- codon()
  #   #File path seems good
  #   df <- read.xlsx(inFile$datapath, colNames = TRUE)
  #   
  #   head(df)
  # })
  # output$NCIS <- renderTable({
  #   req(NCIS())
  #   inFile <- NCIS()
  #   #file path seems good
  #   df <- read.fasta.DNA(inFile$datapath)
  #   head(df)
  # })
  
  #Putting this outside
  
  datasetInput <- reactive({
    req(codon())
    inFile_codon <- codon()
    req(NCIS())
    inFile_NCIS <- NCIS()
    
    # withProgress(message = 'Beginning', detail = "This other thing",
    #              style = "notification", value = NULL, {
    #                Sys.sleep(500) #I think this is seconds!
    #              })
    df <- codonoptimize(
      seq = aa(),
      p = p(),
      n = n(),
      CAI = cai(),
      codons = inFile_codon$datapath,
      nCIS = inFile_NCIS$datapath
    )
    
    df <- df[,c("n", "GC", "CAI","CFD", "nCIS", "nRE", "mrna")]
    
  })
  
  output$CO <- renderTable({
    
    
    datasetInput()
  })
  
  
  # # Downloadable xlsx of selected dataset ----
  # output$downloadData <- downloadHandler(
  #   filename = "sequences.xlsx",
  # 
  #   content = function(file) {
  #     # write.csv(datasetInput(), file, row.names = FALSE)
  # 
  #     write.xlsx(datasetInput(), file, row.names = FALSE)
  #   }
  # )
  # 
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sequences", ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(datasetInput(), file, row.names = FALSE)
    }
  )
  
  

  #Good example
  # output$CO <- renderTable({
  #   req(codon())
  #   inFile_codon <- codon()
  #   req(NCIS())
  #   inFile_NCIS <- NCIS()
  #   
  #   df <- codonoptimize(seq = aa(),
  #                       p = p(),
  #                       n = n(),
  #                       CAI = cai(),
  #                       codons = inFile_codon$datapath,
  #                       nCIS = inFile_NCIS$datapath
  #   )
  #   
  #   df <- df[,c("n", "GC", "CAI","CFD", "nCIS", "nRE", "mrna")]
  # })
  # 
  # # Downloadable xlsx of selected dataset ----
  # output$downloadData <- downloadHandler(
  #   filename = "sequences.xlsx",
  #   
  #   content = function(file) {
  #     # write.csv(datasetInput(), file, row.names = FALSE)
  #     
  #     write.xlsx(df(), file, row.names = FALSE)
  #   }
  # )
  
}

shinyApp(ui, server)