#!/usr/local/bin/Rscript

###
### fista.mat2.weight.S.R
###
### 2017.03.29 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlc.dir = "/home/morii/work/github/moriiism/mchanlc"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mchanlc.dir, "scriptR/fista.mat2.weight.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
infile1     = args[1]
infile2     = args[2]
freq.file   = args[3]
mat.file    = args[4]
mat.w.file  = args[5]
nlambda     = as.integer(args[6])

fista.mat2.weight(infile1, infile2, freq.file, mat.file, mat.w.file, nlambda)
