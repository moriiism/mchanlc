###
### fista.mat2.weight.cv.R
###
### 2017.03.30 M.Morii
### 2017.04.19 M.Morii
###  for 2 light curve
###
### 2017.05.18 M.Morii
###  for cross validation to determine lambda
###  

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlcdir = "/home/morii/work/github/moriiism/mchanlc"

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mitooldir, "script/mirlib/rand.R", sep="/") )
source( paste(mchanlcdir, "scriptR/fista.lib.R", sep="/") )
options(digits=10)
library(grplasso)
library(plotrix)


fista.mat2.weight.cv <- function(infile1, infile2, freq.file, nlambda)
{
    printf(" --- fista.mat ---- \n")
    printf(" infile1    = %s\n", infile1)
    printf(" infile2    = %s\n", infile2)
    printf(" freq.file = %s\n", freq.file)
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
    ## add random index
    ##
    ri.vec1 = GetRandIndex(nrow1)
    ri.vec2 = GetRandIndex(nrow2)
    data1.ri.df = cbind(data1.df, ri.vec1)
    data2.ri.df = cbind(data2.df, ri.vec2)

    ##
    ## Get Nfold Index
    ##

    nelm = nrow1 / nfold
    
    for(iterm in 1:nterm){
        index.st = 



        

        data1.ri.tr.df = data1.ri.df[,]

        
        
        ri.vec1[]

    }



    ##
    ## split into training + test
    ##
    data1.  data1.ri.df
    

    
    ##
    ## make fourier matrix file
    ##

    

}


    

###    ##
###    ## read fourier matrix file
###    ##
###    mat.w.df <- read.table(mat.w.file, comment.char = "#")
###    A.w.mat = data.matrix(mat.w.df)
###    print("matrix filled")
###
###    printf("nrow(A.w.mat) = %d\n", nrow(A.w.mat))
###    printf("ncol(A.w.mat) = %d\n", ncol(A.w.mat))
###
###
###    mat.df <- read.table(mat.file, comment.char = "#")
###    A.mat = data.matrix(mat.df)
###    print("matrix filled")
###
###    B.mat = cbind(1, A.w.mat)
###    h.vec = as.vector(c(data1.df[,2] * sqrt(nrow2), data2.df[,2] * sqrt(nrow1) ))
###
###    ##
###    ## normalize and standardize
###    ##
###    ##mean = mean(h.vec)
###    ##sd   = sd(h.vec)
###    ##h.vec = h.vec - mean
###    ##h.vec = h.vec / sd
###    
###    ##
###    ## calc lambda.array
###    ##
###    index.vec = c(NA, 1:(ncol/4), 1:(ncol/4), 1:(ncol/4), 1:(ncol/4))
###
###    #### adhoc
###    
#####    lambda.array <- lambdamax(B.mat, y = h.vec, index = index.vec,
#####                              penscale = sqrt, model = LinReg()) * 10 * 0.5^(0:(nlambda-2))
###    lambda.array <- lambdamax(B.mat, y = h.vec, index = index.vec,
###                              penscale = sqrt, model = LinReg()) * 200 * 0.8^(0:(nlambda-1))
###    
#####    lambda.array = c(lambda.array, 0.0)
###
###    #### ad hoc
#####    lambda.array = 1e4 * 0.5^(0:(nlambda-2))
#####    lambda.array = c(lambda.array, 0.0)
###    ####
###
###    hist.lo = -0.005
###    hist.up =  0.005
###    hist.nbin = 200
###    hist.delta = (hist.up - hist.lo) / hist.nbin
###    
###    tolerance = 5.e-8
###    n.nonzero = numeric(nlambda)
###    for(ilambda in 1:nlambda){
###        lambda = lambda.array[ilambda]
###        printf("lambda = %e\n", lambda)
###
###        eta = 1.2
###        x = rep(0.0, ncol)
###        x.pre = x
###        y = x
###        L = 1e-3
###        L.pre = L
###        k.max = 500
###        t = 1
###
###        cost = 0.0
###        cost.pre = cost
###        for(k in 1 : k.max){
###            ## printf("k = %d\n", k)
###            ik = FindIk(y, L.pre, eta, A.w.mat, h.vec, lambda)
###            ## printf("ik = %d\n", ik)
###            L = eta**ik * L.pre
###            x = ProxMap(y, L, A.w.mat, h.vec, lambda)
###            t.new = (1. + sqrt(1. + 4 * t**2))/2.
###            y.new = x + (t - 1.)/ t.new * (x - x.pre)
###            x.pre = x
###            y = y.new
###            t = t.new
###            L.pre = L
###
###            cost = FGFunc(y.new, A.w.mat, h.vec, lambda)
###            ##printf("ik = %d, L = %e, t = %e, FGFunc = %e\n",
###            ##       ik, L, t, FGFunc(y.new, A.w.mat, h.vec, lambda))
###
###            ##printf("diff = %e\n", (cost.pre - cost) / cost)
###            
###            if(k > 1 && (cost.pre - cost) / cost < tolerance){
###                printf("k = %d, cost = %e\n", k, cost)
###                break
###            }
###            cost.pre = cost
###        }
###
###        n.nonzero[ilambda] = 0
###        for(i in 1:length(x)){
###            if(abs(x[i]) > lambda){
###                n.nonzero[ilambda] = n.nonzero[ilambda] + 1
###            }
###        }
###        printf("n.nonzero = %d\n", n.nonzero[ilambda])
###
###        ##
###        ## delta.time
###        ##
###        N = ncol/4
###        delta.time12.vec = numeric(0)
###        weight.delta.time12.vec = numeric(0)
###
###        delta.time12.vec.tmp = numeric(0)
###        weight.delta.time12.vec.tmp = numeric(0)
###        period.vec.tmp = numeric(0)
###        
###        delta.time1.vec = numeric(0)
###        delta.time2.vec = numeric(0)
###        isel.vec = numeric(0)
###
###        x.common = rep(0.0, ncol)
###        
###        for(i in 1:N){
###            if(abs(x[i]) > 0.0  && abs(x[i + N]) > 0.0 && abs(x[i + 2 * N]) > 0.0 && abs(x[i + 3 * N]) > 0.0 ) {
###                freq = freq.lo + delta.freq * (i - 0.5)
###                omega = 2 * pi * freq
###                delta.time1 = GetAlphaOfSinCos(x[i + N], x[i]) / omega
###                delta.time2 = GetAlphaOfSinCos(x[i + 3 * N], x[i + 2 * N]) / omega
###                delta.time1 = GetDeltaTime(delta.time1, 1./freq)
###                delta.time2 = GetDeltaTime(delta.time2, 1./freq)
###                delta.time12 = delta.time1 - delta.time2
###                weight.delta.time12 = sqrt( x[i]**2 + x[i + N]**2 + x[i + 2 * N]**2 + x[i + 3 * N]**2 )
###
###                ## select common between X and Opt
###                if(0.0 < delta.time12 && delta.time12 < 0.01){
###                    x.common[i] = x[i]
###                    x.common[i + N] = x[i + N]
###                    x.common[i + 2 * N] = x[i + 2 * N]
###                    x.common[i + 3 * N] = x[i + 3 * N]
###                }
###                
###                ##printf("delta.time1 = %e\n", delta.time1)
###                ##printf("delta.time2 = %e\n", delta.time2)                
###                ##printf("delta.time12 = %e\n", delta.time12)
###                ##printf("T = %e\n", 1./freq)
###                ##printf("weight.delta.time12 = %e\n", weight.delta.time12)
###
###                delta.time12.vec.tmp = c(delta.time12.vec.tmp, delta.time1 - delta.time2)
###                weight.delta.time12.vec.tmp = c(weight.delta.time12.vec.tmp, weight.delta.time12)
###                period.vec.tmp = c(period.vec.tmp, 1./freq)
###                
###                if(hist.lo < delta.time12 && delta.time12 < hist.up){
###                    delta.time12.vec = c(delta.time12.vec, delta.time1 - delta.time2)
###                    weight.delta.time12.vec = c(weight.delta.time12.vec, weight.delta.time12)
###                    ##weight.delta.time12.vec = c(weight.delta.time12.vec, 1.0)
###                    
###                    ## ibin = GetIbin(delta.time1 - delta.time2, hist.lo, hist.up, hist.nbin)
###                    ## hist.oval[ibin] = hist.oval[ibin] + 1.0
###                    
###                }
###                isel.vec = c(isel.vec, i)
###            }
###        }
###        hist(delta.time12.vec, breaks=seq(hist.lo,hist.up,hist.delta))
###        ## weighted.hist(delta.time12.vec, weight.delta.time12.vec, breaks=seq(hist.lo,hist.up,hist.delta))
###
###
####        ## adhoc
####        for(isel in 1:length(isel.vec)){
####            xsel = rep(0.0, length(x))
####            N = length(x)/4
####            index = isel.vec[isel]
####            xsel[index] = x[index]
####            xsel[index + N] = x[index + N]
####            xsel[index + 2 * N] = x[index + 2 * N]
####            xsel[index + 3 * N] = x[index + 3 * N]
####            h.rec.sel.vec = A.w.mat %*% xsel
####            lc.sel.rec = cbind(data.df[, 1], h.rec.sel.vec)
####            outfile = sprintf("rec_%2.2d_sel_%2.2d.dat", ilambda, isel)
####            write(t(lc.sel.rec), file=outfile, ncolumns = 2)
####        }
###
###        ##
###        ## lag v.s. norm
###        ##
###        lag.norm = cbind(delta.time12.vec.tmp, weight.delta.time12.vec.tmp)
###        outfile = sprintf("lag_norm_%2.2d.dat", ilambda)
###        write(t(lag.norm), file=outfile, ncolumns = 2)
###
###        ##
###        ## lag v.s. norm v.s. period
###        ##
###        lag.norm = cbind(delta.time12.vec.tmp, weight.delta.time12.vec.tmp, period.vec.tmp)
###        outfile = sprintf("lag_norm_period_%2.2d.dat", ilambda)
###        write(t(lag.norm), file=outfile, ncolumns = 3)
###
###        
###       
###        h.rec.vec = A.mat %*% x
###
###        ##
###        ## h.rec.vec = h.rec.vec * sd + mean
###        ##
###        
###        lc.rec = cbind(data.df[, 1], h.rec.vec)
###        outfile = sprintf("rec_%2.2d.dat", ilambda)
###        write(t(lc.rec), file=outfile, ncolumns = 2)
###
###        outfile = sprintf("x_%2.2d.dat", ilambda)
###        x.mat = cbind(c(1:(length(x)/2),1:(length(x)/2)), x)
###        write(t(x.mat), file=outfile, ncolumns = 2)
###
###        ##
###        ## common
###        ##
###        h.rec.common.vec = A.mat %*% x.common
###
###        lc.common.rec = cbind(data.df[, 1], h.rec.common.vec)
###        outfile = sprintf("rec_common.%2.2d.dat", ilambda)
###        write(t(lc.common.rec), file=outfile, ncolumns = 2)
###
###        outfile = sprintf("x_common.%2.2d.dat", ilambda)
###        x.common.mat = cbind(c(1:(length(x.common)/2),1:(length(x.common)/2)), x.common)
###        write(t(x.common.mat), file=outfile, ncolumns = 2)
###
###        power.common = 0.0
###        for(i in 1:N){
###            power.common = power.common + x.common[i]**2 + x.common[i + N]**2 + x.common[i + 2 * N]**2 + x.common[i + 3 * N]**2
###        }
###        printf("power.common = %e\n", power.common)
###
###        power = 0.0
###        for(i in 1:N){
###            power = power + x[i]**2 + x[i + N]**2 + x[i + 2 * N]**2 + x[i + 3 * N]**2
###        }
###        printf("power = %e\n", power)
###
###        
###        ##
###        ## uncommon 
###        ##
###                
###        lc.uncommon.rec = cbind(data.df[, 1], h.rec.vec - h.rec.common.vec)
###        outfile = sprintf("rec_uncommon.%2.2d.dat", ilambda)
###        write(t(lc.uncommon.rec), file=outfile, ncolumns = 2)
###        
###    }
###
###    for(ilambda in 1:nlambda){
###        lambda = lambda.array[ilambda]
###        printf("%d: lambda = %e, n.nonzero[ilambda] = %d \n",
###               ilambda, lambda, n.nonzero[ilambda])
###    }
###    
###}
###
###
###
