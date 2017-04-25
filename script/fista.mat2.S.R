#!/usr/local/bin/Rscript

###
### fista.mat2.S.R
###
### 2017.03.29 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlc.dir = "/home/morii/work/github/moriiism/mchanlc"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mchanlc.dir, "scriptR/fista.mat2.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
infile1     = args[1]
infile2     = args[2]
freq.file  = args[3]
mat.file   = args[4]
nlambda    = as.integer(args[5])

fista.mat2(infile1, infile2, freq.file, mat.file, nlambda)
