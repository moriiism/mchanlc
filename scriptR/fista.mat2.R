###
### fista.mat2.R
###
### 2017.03.30 M.Morii
### 2017.04.19 M.Morii
###  for 2 light curve
###

mitooldir = "/home/morii/work/github/moriiism/mitool"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
options(digits=10)
library(grplasso)
library(plotrix)

###
### a sin(theta) + b cos(theta) = sqrt(a**2 + b**2) * sin(theta + alpha)
###
GetAlphaOfSinCos <- function(a, b)
{
    dist = sqrt(a**2 + b**2)
    cos.alpha = a / dist
    sin.alpha = b / dist
    if(sin.alpha >= 0.0){
        alpha = acos(cos.alpha)
    }
    else {
        alpha = - acos(cos.alpha)
    }
    return(alpha)
}

###
### 0 <= delta.time <= period
###
GetDeltaTime <- function(time, period)
{
    delta.time = time %% period
}




GetIbin <- function(xval, xlo, xup, nbinx)
{
    delta.xval = (xup - xlo) / nbinx
    ibin = ceiling((xval - xlo) / delta.xval)
    return(ibin)
}

###SoftThres <- function(x, lambda)
###{
###    if(lambda < x){
###        ans = x - lambda
###    } else if(-lambda <= x && x <= lambda){
###        ans = 0.0
###    } else {
###        ans = x + lambda
###    }
###    return(ans)
###}
###
###SoftThres2 <- function(a, b, L, lambda)
###{
###    len = sqrt(a**2 + b**2)
###    ans = 1./L / len * SoftThres(L * len, lambda)
###    return(ans)
###}

SoftThres2 <- function(a, b, c, d, L, lambda)
{
    l = sqrt(a**2 + b**2 + c**2 + d**2) * L
    ans = 0.0
    if(lambda < l){
        ans = 1.0 - lambda / l
    }
    if(-lambda <= l && l <= lambda){
        ans = 0.0
    }
    if(l < -lambda){
        ans = 1.0 + lambda / l
    }
    return(ans)
}


ProxMap <- function(y, L, A.mat, h.vec, lambda)
{
    b = y - 2.0/L * t(A.mat) %*% (A.mat %*% y - h.vec)
    N2 = length(b)
    N = N2 / 2
    b.new = rep(0.0, 2 * N)
    for(i in 1:(N/2)){
        factor = SoftThres2(b[i], b[i + N/2], b[i + N], b[i + N/2 * 3], L, lambda)
        b.new[i]       = b[i] * factor
        b.new[i + N/2] = b[i + N/2] * factor
        b.new[i + N]   = b[i + N] * factor
        b.new[i + N/2 * 3] = b[i + N/2 * 3] * factor
    }
    return(b.new)
}

GFunc <- function(J.vec, lambda){
    N2 = length(J.vec)
    N = N2 / 2
    sum = 0.0
    for(i in 1:(N/2)){
        sum = sum + sqrt(J.vec[i]**2 + J.vec[i + N/2]**2 + J.vec[i + N]**2 + J.vec[i + N/2 * 3]**2)
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
    ik.max = 1000
    ik = 0
    while(ik <= ik.max){
        L = eta**ik * L.pre
        pLy = ProxMap(y, L, A.mat, h.vec, lambda)
        FGFunc.val = FGFunc(pLy, A.mat, h.vec, lambda)
        QFunc.val  = QFunc(pLy, y, L, A.mat, h.vec, lambda)
        ### printf("FGFunc.val = %e, QFunc.val = %e\n", FGFunc.val, QFunc.val)
        if(FGFunc.val <= QFunc.val){
            break
        }
        ik = ik + 1
    }
    return(ik)
}

fista.mat2 <- function(infile1, infile2, freq.file, mat.file, nlambda)
{
    printf(" --- fista.mat ---- \n")
    printf(" infile1    = %s\n", infile1)
    printf(" infile2    = %s\n", infile2)
    printf(" freq.file = %s\n", freq.file)
    printf(" mat.file  = %s\n", mat.file)
    printf(" nlambda   = %d\n", nlambda)
    printf(" --- fista.mat ---- \n")

    ##
    ## read infile1
    ##
    data1.df <- read.table(infile1, comment.char = "#")
    printf("ncol(data1.df) = %d\n", ncol(data1.df))
    printf("nrow(data1.df) = %d\n", nrow(data1.df))
    nrow1 = nrow(data1.df)

    ##
    ## read infile2
    ##
    data2.df <- read.table(infile2, comment.char = "#")
    printf("ncol(data2.df) = %d\n", ncol(data2.df))
    printf("nrow(data2.df) = %d\n", nrow(data2.df))
    nrow2 = nrow(data2.df)

    data.df = rbind(data1.df, data2.df)

    ##
    ## read frequency setup file
    ##
    freq.df <- read.table(freq.file, comment.char = "#")
    freq.lo = freq.df[1, 1]
    freq.up = freq.df[1, 2]
    nfreq   = freq.df[1, 3]

    delta.freq = (freq.up - freq.lo) / nfreq
    ncol = nfreq * 2 * 2

    ##
    ## read fourier matrix file
    ##
    mat.df <- read.table(mat.file, comment.char = "#")
    A.mat = data.matrix(mat.df)
    print("matrix filled")

    printf("nrow(A.mat) = %d\n", nrow(A.mat))
    printf("ncol(A.mat) = %d\n", ncol(A.mat))


    B.mat = cbind(1, A.mat)
    h.vec = as.vector(data.df[,2])

    ##
    ## normalize and standardize
    ##
    ##mean = mean(h.vec)
    ##sd   = sd(h.vec)
    ##h.vec = h.vec - mean
    ##h.vec = h.vec / sd
    
    ##
    ## calc lambda.array
    ##
    index.vec = c(NA, 1:(ncol/4), 1:(ncol/4), 1:(ncol/4), 1:(ncol/4))

    #### adhoc
    
##    lambda.array <- lambdamax(B.mat, y = h.vec, index = index.vec,
##                              penscale = sqrt, model = LinReg()) * 10 * 0.5^(0:(nlambda-2))
    lambda.array <- lambdamax(B.mat, y = h.vec, index = index.vec,
                              penscale = sqrt, model = LinReg()) * 10 * 0.5^(0:(nlambda-1))
    
##    lambda.array = c(lambda.array, 0.0)

    #### ad hoc
##    lambda.array = 1e4 * 0.5^(0:(nlambda-2))
##    lambda.array = c(lambda.array, 0.0)
    ####

    hist.lo = -0.005
    hist.up =  0.005
    hist.nbin = 200
    hist.delta = (hist.up - hist.lo) / hist.nbin
    
    tolerance = 5.e-8
    n.nonzero = numeric(nlambda)
    for(ilambda in 1:nlambda){
        lambda = lambda.array[ilambda]
        printf("lambda = %e\n", lambda)

        eta = 1.2
        x = rep(0.0, ncol)
        x.pre = x
        y = x
        L = 1e-3
        L.pre = L
        k.max = 500
        t = 1

        cost = 0.0
        cost.pre = cost
        for(k in 1 : k.max){
            ## printf("k = %d\n", k)
            ik = FindIk(y, L.pre, eta, A.mat, h.vec, lambda)
            ## printf("ik = %d\n", ik)
            L = eta**ik * L.pre
            x = ProxMap(y, L, A.mat, h.vec, lambda)
            t.new = (1. + sqrt(1. + 4 * t**2))/2.
            y.new = x + (t - 1.)/ t.new * (x - x.pre)
            x.pre = x
            y = y.new
            t = t.new
            L.pre = L

            cost = FGFunc(y.new, A.mat, h.vec, lambda)
            ##printf("ik = %d, L = %e, t = %e, FGFunc = %e\n",
            ##       ik, L, t, FGFunc(y.new, A.mat, h.vec, lambda))

            ##printf("diff = %e\n", (cost.pre - cost) / cost)
            
            if(k > 1 && (cost.pre - cost) / cost < tolerance){
                printf("k = %d, cost = %e\n", k, cost)
                break
            }
            cost.pre = cost
        }

        n.nonzero[ilambda] = 0
        for(i in 1:length(x)){
            if(abs(x[i]) > lambda){
                n.nonzero[ilambda] = n.nonzero[ilambda] + 1
            }
        }
        printf("n.nonzero = %d\n", n.nonzero[ilambda])

        ##
        ## delta.time
        ##
        N = ncol/4
        delta.time12.vec = numeric(0)
        weight.delta.time12.vec = numeric(0)
        delta.time1.vec = numeric(0)
        delta.time2.vec = numeric(0)
        isel.vec = numeric(0)
        for(i in 1:N){
            if(abs(x[i]) > 0.0  && abs(x[i + N]) > 0.0 && abs(x[i + 2 * N]) > 0.0 && abs(x[i + 3 * N]) > 0.0 ) {
                freq = freq.lo + delta.freq * (i - 0.5)
                omega = 2 * pi * freq
                delta.time1 = GetAlphaOfSinCos(x[i + N], x[i]) / omega
                delta.time2 = GetAlphaOfSinCos(x[i + 3 * N], x[i + 2 * N]) / omega
                delta.time1 = GetDeltaTime(delta.time1, 1./freq)
                delta.time2 = GetDeltaTime(delta.time2, 1./freq)
                delta.time12 = delta.time1 - delta.time2
                weight.delta.time12 = sqrt( x[i]**2 + x[i + N]**2 + x[i + 2 * N]**2 + x[i + 3 * N]**2 )

                printf("delta.time1 = %e\n", delta.time1)
                printf("delta.time2 = %e\n", delta.time2)                
                printf("delta.time12 = %e\n", delta.time12)
                printf("T = %e\n", 1./freq)

                printf("weight.delta.time12 = %e\n", weight.delta.time12)
                
                if(hist.lo < delta.time12 && delta.time12 < hist.up){
                    delta.time12.vec = c(delta.time12.vec, delta.time1 - delta.time2)
                    weight.delta.time12.vec = c(weight.delta.time12.vec, weight.delta.time12)
                    ##weight.delta.time12.vec = c(weight.delta.time12.vec, 1.0)
                    
                    ## ibin = GetIbin(delta.time1 - delta.time2, hist.lo, hist.up, hist.nbin)
                    ## hist.oval[ibin] = hist.oval[ibin] + 1.0
                    
                }
                isel.vec = c(isel.vec, i)
            }
        }
        hist(delta.time12.vec, breaks=seq(hist.lo,hist.up,hist.delta))
        ## weighted.hist(delta.time12.vec, weight.delta.time12.vec, breaks=seq(hist.lo,hist.up,hist.delta))


#        ## adhoc
#        for(isel in 1:length(isel.vec)){
#            xsel = rep(0.0, length(x))
#            N = length(x)/4
#            index = isel.vec[isel]
#            xsel[index] = x[index]
#            xsel[index + N] = x[index + N]
#            xsel[index + 2 * N] = x[index + 2 * N]
#            xsel[index + 3 * N] = x[index + 3 * N]
#            h.rec.sel.vec = A.mat %*% xsel
#            lc.sel.rec = cbind(data.df[, 1], h.rec.sel.vec)
#            outfile = sprintf("rec_%2.2d_sel_%2.2d.dat", ilambda, isel)
#            write(t(lc.sel.rec), file=outfile, ncolumns = 2)
#        }
        
        h.rec.vec = A.mat %*% x

        ##
        ## h.rec.vec = h.rec.vec * sd + mean
        ##
        
        lc.rec = cbind(data.df[, 1], h.rec.vec)
        outfile = sprintf("rec_%2.2d.dat", ilambda)
        write(t(lc.rec), file=outfile, ncolumns = 2)

        outfile = sprintf("x_%2.2d.dat", ilambda)
        x.mat = cbind(c(1:(length(x)/2),1:(length(x)/2)), x)
        write(t(x.mat), file=outfile, ncolumns = 2)
    }

    for(ilambda in 1:nlambda){
        lambda = lambda.array[ilambda]
        printf("%d: lambda = %e, n.nonzero[ilambda] = %d \n",
               ilambda, lambda, n.nonzero[ilambda])
    }
    
}



