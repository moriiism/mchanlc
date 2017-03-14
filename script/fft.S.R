#!/usr/local/bin/Rscript

###
### fft.S.R
###
### 2017.03.13 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )

###

args <- commandArgs(TRUE)
infile   = args[1]

data.df <- read.table(infile, comment.char = "#")
printf("ncol(data.df) = %d\n", ncol(data.df))
printf("nrow(data.df) = %d\n", nrow(data.df))

time.lo = data.df[1,]
    



