#!/usr/local/bin/Rscript

###
### sim.trig.S.R
###
### Simulate light curve of linear combination of
### trigonometric functions (sine function),
### using frequency and phase information, on GTI duration.
###
### 2017.03.31 M.Morii
###
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlc.dir = "/home/morii/work/github/moriiism/mchanlc"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mchanlc.dir, "scriptR/sim.trig.R", sep="/") )
options(digits=10)

###
### main
###

args <- commandArgs(TRUE)
gti.info.file  = args[1]
trig.info.file = args[2]
delta.time     = as.numeric(args[3])
outfile        = args[4]

###
### read gti info file
###
gti.df <- read.table(gti.info.file, comment.char = "#")
ngti = ncol(gti.df)

print(gti.df)

###
### read trig info file
###
trig.df <- read.table(trig.info.file, comment.char = "#")

print(trig.df)

###
### make gti.itime
###

time.lo = gti.df[1, 1]
gti.itime.df = GtiItime(time.lo, delta.time, gti.df)

print(gti.itime.df)

time.vec = numeric(0)
oval.vec = numeric(0)
ivec = 0
for(igti in 1:nrow(gti.itime.df)){
    itime.st = gti.itime.df[igti, 1]
    itime.ed = gti.itime.df[igti, 2]
    for(itime in itime.st:itime.ed){
        time = time.lo + (itime + 0.5) * delta.time
        time.vec[ivec] = time
        oval.vec[ivec] = TrigFunc(time, trig.df)
        ivec = ivec + 1
    }
}

lc = cbind(time.vec, oval.vec)
write(t(lc), file=outfile, ncolumns = 2)

