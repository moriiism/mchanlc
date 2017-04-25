#!/usr/local/bin/Rscript

###
### grplasso.each.S.R
###
### 2017.03.30 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlc.dir = "/home/morii/work/github/moriiism/mchanlc"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mchanlc.dir, "scriptR/grplasso.each.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
infile     = args[1]
freq.file  = args[2]
mat.file   = args[3]
nlambda    = as.integer(args[4])

GrpLasso(infile, freq.file, mat.file, nlambda)

