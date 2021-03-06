% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_write_fasta.R
\name{read_fasta}
\alias{read_fasta}
\title{Read sequences in FASTA format}
\usage{
read_fasta(file, type = "DNA")
}
\arguments{
\item{file}{A path to a file, a connection or a character string.}

\item{type}{Type of data. Can be "DNA" (the default), "RNA" or "AA".}
}
\value{
A DNA, RNA or AA vector (depending on \code{type} argument).
}
\description{
Read sequences in FASTA format
}
\seealso{
Other input/output operations: 
\code{\link{write_fasta}()}
}
\concept{input/output operations}
