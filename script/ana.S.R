#!/usr/local/bin/Rscript

###
### ana.S.R
###
### 2017.02.15 M.Morii
###
###

library(grplasso)

args <- commandArgs(TRUE)

infile   = args[1]

data.df <- read.table(infile, comment.char = "#")
print(ncol(data.df))
print(nrow(data.df))
N = nrow(data.df)

freq.lo =  5.e-2
freq.up =  9.e+4
nfreq = 5000
delta.freq = (freq.up - freq.lo) / nfreq
A.mat = matrix(0.0, nrow=N, ncol=nfreq)

for(i in 1:(N/2)){
    for(n in 1:(nfreq/2)){
        freq = freq.lo + delta.freq * (n - 0.5)
        A.mat[i, n] = 2 * delta.freq * cos(2 * pi * freq * data.df[i, 1])
        A.mat[N/2 + i, nfreq/2 + n] = 2 * delta.freq * sin(2 * pi * freq * data.df[i, 1])
    }
}
print("matrix filled")

B.mat = cbind(1, A.mat)
h.vec = as.vector(data.df[,2])
index.vec = c(NA, 1:(nfreq/2), 1:(nfreq/2))

lambda <- lambdamax(B.mat, y = h.vec, index = index.vec,
                    penscale = sqrt, model = LinReg()) * 0.5^(0:10)

# print(head(B.mat))
# print(head(h.vec))
# print(head(index.vec))


fit = grplasso(x = B.mat, y = h.vec, index = index.vec, model = LinReg(),
    lambda = lambda, penscale = sqrt, center = TRUE,
    standardize = TRUE,
    control = grpl.control(update.hess = "lambda", trace = 0))

print(fit)

# coef(fit)

png("plot.png")
plot(fit)
dev.off()

write(coef(fit), file="fit.dat")


# 57190.553103 12.543 
# 57208.453930 14.890


res.vec = coef(fit)[,10]

ncol(A.mat)
nrow(A.mat)
length(res.vec)

h.rec.vec = B.mat %*% res.vec
lc.rec = cbind(data.df[, 1], h.rec.vec)
write(t(lc.rec), file="rec.dat", ncolumns = 2)
