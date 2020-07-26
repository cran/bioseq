## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bioseq)

## -----------------------------------------------------------------------------
x <- dna(Seq_1 = "ACCTAG", Seq_2 = "GGTATATACC", Seq_3 = "AGTC")
is_dna(x)
x

## -----------------------------------------------------------------------------
x[c("Seq_3", "Seq_1")]
x[2]
x[c(FALSE, FALSE, TRUE)]

## -----------------------------------------------------------------------------
y <- dna("?AcGF")
y

## -----------------------------------------------------------------------------
x_dna <- dna("ATGTCACCACAAACAGAGACT")
x_dna

x_rna <- seq_transcribe(x_dna)
x_rna

x_aa <- seq_translate(x_rna)
x_aa

## -----------------------------------------------------------------------------
dna_from_rna <- seq_rev_transcribe(x_rna)
dna_from_rna

dna_from_aa <- seq_rev_translate(x_aa)
dna_from_aa

## -----------------------------------------------------------------------------
x_dna_comp <- seq_complement(x_dna)
x_dna_comp_rev <- seq_reverse(x_dna_comp)

dna(x_dna, x_dna_comp, x_dna_comp_rev)

## -----------------------------------------------------------------------------
x <- dna("CTGAAAACTG", "ATGAAAACTG", "CTGCTG")

## -----------------------------------------------------------------------------
x[seq_detect_pattern(x, "AAAA")]

## -----------------------------------------------------------------------------
x[seq_detect_pattern(x, "A{4}")]

## -----------------------------------------------------------------------------
# This works
x[seq_detect_pattern(x, dna("AAAA"))]

## ---- eval=FALSE--------------------------------------------------------------
#  # This fails because x is a DNA vector and pattern is an amino acid vector
#  x[seq_detect_pattern(x, aa("AAAA"))]

## -----------------------------------------------------------------------------
# This works because W can be A or T.
x[seq_detect_pattern(x, dna("WAWA"))]

## -----------------------------------------------------------------------------
seq_disambiguate_IUPAC(dna("WAWA"))

## -----------------------------------------------------------------------------
seq_remove_pattern(x, "A{4}")

## -----------------------------------------------------------------------------
seq_replace_pattern(x,
                    pattern = dna("AAAA"),
                    replacement = dna("----"))

## -----------------------------------------------------------------------------
x <- seq_remove_pattern(x, "A{4}")
seq_replace_position(x, 4, 6,
                     replacement = dna("CCC"))

## -----------------------------------------------------------------------------
x
seq_replace_position(x, 1:3, 6,
                     replacement = dna("-", "--", "---"))

