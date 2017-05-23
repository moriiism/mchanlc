#!/usr/local/bin/Rscript

###
### fista.mat2.weight.slide.S.R
###
### 2017.05.23 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlc.dir = "/home/morii/work/github/moriiism/mchanlc"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mchanlc.dir, "scriptR/fista.mat2.weight.slide.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
infile1     = args[1]
infile2     = args[2]
freq.file   = args[3]
lambda      = as.numeric(args[4])
slide.file  = args[5]
outdir      = args[6]

fista.mat2.weight.slide(infile1, infile2, freq.file, lambda, slide.file, outdir)
