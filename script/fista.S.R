#!/usr/local/bin/Rscript

###
### fista.S.R
###
### 2017.03.13 M.Morii
###
###

options(digits=10)
mitooldir = "/home/morii/work/github/moriiism/mitool"

library(grplasso)

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
### library(grplasso)

###
SoftThres <- function(x, lambda)
{
    if(lambda < x){
        ans = x - lambda
    } else if(-lambda <= x && x <= lambda){
        ans = 0.0
    } else {
        ans = x + lambda
    }
    return(ans)
}

SoftThres2 <- function(a, b, L, lambda)
{
    len = sqrt(a**2 + b**2)
    ans = 1./L / len * SoftThres(L * len, lambda)
    return(ans)
}

ProxMap <- function(y, L, A.mat, h.vec, lambda)
{
    b = y - 2.0/L * t(A.mat) %*% (A.mat %*% y - h.vec)
    N = length(b)
    b.new = rep(0.0, N)
    for(i in 1:(N/2)){
        factor = SoftThres2(b[i], b[i + N/2], L, lambda)
        b.new[i] = b[i] * factor
        b.new[i + N/2] = b[i + N/2] * factor
    }
    return(b.new)
}

GFunc <- function(J.vec, lambda){
    N = length(J.vec)
    sum = 0.0    
    for(i in 1:(N/2)){
        sum = sum + sqrt(J.vec[i]**2 + J.vec[i + N/2]**2)
    }
    ans = lambda * sum
    return(ans)
}

DiffFFunc <- function(y, A.mat, h.vec){
    ans = 2 * t(A.mat) %*% (A.mat %*% y - h.vec)
    return(ans)
}

FFunc <- function(J.vec, A.mat, h.vec){
    vec.tmp = h.vec - A.mat %*% J.vec
    ans = sum( vec.tmp * vec.tmp )
    return(ans)
}

QFunc <- function(x, y, L, A.mat, h.vec, lambda){
    term1 = FFunc(y, A.mat, h.vec)
    term2 = sum((x - y) * DiffFFunc(y, A.mat, h.vec))
    term3 = L / 2. * sum((x - y) * (x - y))
    term4 = GFunc(x, lambda)
    ans = term1 + term2 + term3 + term4
    return(ans)
}

FGFunc <- function(J.vec, A.mat, h.vec, lambda){
    ans = FFunc(J.vec, A.mat, h.vec) + GFunc(J.vec, lambda)
    return(ans)
}

FindIk <- function(y, L.pre, eta, A.mat, h.vec, lambda){
    ik.max = 100
    ik = 1
    for(ik in 1 : ik.max){
        L = eta**ik * L.pre
        pLy = ProxMap(y, L, A.mat, h.vec, lambda)
        FGFunc.val = FGFunc(pLy, A.mat, h.vec, lambda)
        QFunc.val = QFunc(pLy, y, L, A.mat, h.vec, lambda)
        if(FGFunc.val <= QFunc.val){
            break
        }
    }
    printf("ik = %d\n", ik)
    
    return(ik)
}

###

args <- commandArgs(TRUE)
infile   = args[1]

data.df <- read.table(infile, comment.char = "#")
printf("ncol(data.df) = %d\n", ncol(data.df))
printf("nrow(data.df) = %d\n", nrow(data.df))
nrow = nrow(data.df)


freq.lo =  0.1
# freq.up =  9.e+4
freq.up =  1.0
nfreq = 100
delta.freq = (freq.up - freq.lo) / nfreq
ncol = nfreq * 2
A.mat = matrix(0.0, nrow=nrow, ncol=ncol)

for(irow in 1:nrow){
    for(icol in 1:(ncol/2)){
        freq = freq.lo + delta.freq * (icol - 0.5)
        A.mat[irow, icol]          = 2 * delta.freq * cos(2 * pi * freq * data.df[irow, 1])
        A.mat[irow, ncol/2 + icol] = 2 * delta.freq * sin(2 * pi * freq * data.df[irow, 1])
    }
}
print("matrix filled")

### calc lambda.array
B.mat = cbind(1, A.mat)
h.vec = as.vector(data.df[,2])
index.vec = c(NA, 1:(ncol/2), 1:(ncol/2))

lambda <- lambdamax(B.mat, y = h.vec, index = index.vec,
                    penscale = sqrt, model = LinReg()) * 0.5^(0:10)

lambda


eta = 1.2
h.vec = data.df[1:nrow,2]
x = rep(0.0, ncol)
x.pre = x
y = x
L = 1e-5
L.pre = L
k.max = 10
t = 1
lambda = 4.862804651

for(k in 1 : k.max){
    printf("k = %d\n", k)
    ik = FindIk(y, L.pre, eta, A.mat, h.vec, lambda)
    L = eta**ik * L.pre
    x = ProxMap(y, L, A.mat, h.vec, lambda)
    t.new = (1. + sqrt(1. + 4 * t**2))/2.
    y.new = x + (t - 1.)/ t.new * (x - x.pre)


    printf("FGFunc = %e\n", FGFunc(y.new, A.mat, h.vec, lambda))

    x.pre = x
    y = y.new
    t = t.new
    L.pre = L
    printf("L = %e\n", L)
    printf("t = %e\n", t)
###    print(x)
    
}

h.rec.vec = A.mat %*% x
lc.rec = cbind(data.df[, 1], h.rec.vec)
write(t(lc.rec), file="rec.dat", ncolumns = 2)



