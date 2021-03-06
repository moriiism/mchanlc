#!/usr/local/bin/Rscript

###
### fista.mat.S.R
###
### 2017.03.29 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlc.dir = "/home/morii/work/github/moriiism/mchanlc"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mchanlc.dir, "scriptR/fista.mat.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
infile     = args[1]
freq.file  = args[2]
mat.file   = args[3]
nlambda    = as.integer(args[4])

fista.mat(infile, freq.file, mat.file, nlambda)
