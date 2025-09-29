## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5
)

## -----------------------------------------------------------------------------
library(bioseq)

## ----eval=FALSE---------------------------------------------------------------
# library(tidyverse)

## ----include=FALSE------------------------------------------------------------
library(tibble)
library(dplyr)
library(stringr)

## -----------------------------------------------------------------------------
data(fragilaria, package = "bioseq")
fra_data <- read_fasta(fragilaria)

fra_data

## -----------------------------------------------------------------------------
seq_nchar(fra_data) %>% range()

## -----------------------------------------------------------------------------
fra_data <- tibble(label = names(fra_data), sequence = fra_data)

## -----------------------------------------------------------------------------
fra_data <- fra_data %>% 
  mutate(genbank_id = str_extract(label, "([^\\s]+)"),
         taxa = str_extract(label, "(?<= ).*")) %>% 
  select(genbank_id, taxa, sequence)

fra_data <- fra_data %>% 
  mutate(n_base = seq_nchar(sequence))

fra_data

## -----------------------------------------------------------------------------
FWD <- dna("AGGTGAAGTAAAAGGTTCWTACTTAAA",
           "AGGTGAAGTTAAAGGTTCWTAYTTAAA",
           "AGGTGAAACTAAAGGTTCWTACTTAAA")

REV <- dna("CAGTWGTWGGTAAATTAGAAGG",
           "CTGTTGTWGGTAAATTAGAAGG")

## -----------------------------------------------------------------------------
seq_disambiguate_IUPAC(FWD)

## -----------------------------------------------------------------------------
fra_data <- fra_data %>% 
  mutate(barcode = seq_crop_pattern(sequence,
                                    pattern_in = list(FWD),
                                    pattern_out = list(REV)))

fra_data

## -----------------------------------------------------------------------------
fra_data <- fra_data %>% 
  filter(seq_nchar(barcode) == 312)

fra_data

## -----------------------------------------------------------------------------
fra_consensus <- fra_data %>% 
  group_by(taxa) %>% 
  summarise(consensus_barcode = seq_consensus(barcode))

fra_consensus

## -----------------------------------------------------------------------------
fra_consensus %>% 
  as_DNAbin(consensus_barcode, taxa) %>% 
  ape::dist.dna() %>% 
  ape::bionj() %>% 
  plot()

## -----------------------------------------------------------------------------
duplicated(fra_consensus$consensus_barcode)

## -----------------------------------------------------------------------------
fra_consensus <- 
  fra_consensus %>% 
  mutate(cluster = seq_cluster(consensus_barcode,
                               threshold = 0.001))
fra_consensus

## -----------------------------------------------------------------------------
fra_consensus <-
  fra_consensus %>% 
  group_by(cluster) %>% 
  summarise(taxa_group = paste(taxa, collapse = "/"),
            consensus_barcode = seq_consensus(consensus_barcode))

fra_consensus %>% 
  as_DNAbin(consensus_barcode, taxa_group) %>% 
  ape::dist.dna() %>% 
  ape::bionj() %>% 
  plot()

## ----eval=FALSE---------------------------------------------------------------
# fra_consensus %>%
#   select(taxa_group, consensus_barcode) %>%
#   deframe() %>%
#   write_fasta("my_sequences.fasta")

