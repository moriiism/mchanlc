#!/usr/local/bin/Rscript

###
### sim.S.R
###
### 2017.03.29 M.Morii
###
###

options(digits=10) 

args <- commandArgs(TRUE)

infile   = args[1]
outfile  = args[2]

data.df <- read.table(infile, comment.char = "#")
print(ncol(data.df))
print(nrow(data.df))
N = nrow(data.df)

freq.1  = 0.1
freq.2  = 0.2
freq.3  = 0.3
norm.1  = 1.0
norm.2  = 2.0
norm.3  = 3.0

h.vec = numeric(N)

for(i in 1:N){
    h.vec[i] =
        norm.1 * sin(2 * pi * freq.1 * data.df[i, 1]) +
            norm.2 * cos(2 * pi * freq.2 * data.df[i, 1]) +
                norm.3 * cos(2 * pi * freq.3 * data.df[i, 1])
}

lc = cbind(data.df[, 1], h.vec)
write(t(lc), file=outfile, ncolumns = 2)
