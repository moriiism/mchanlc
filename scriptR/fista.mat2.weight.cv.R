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
source( paste(mitooldir, "script/mirlib/binning.R", sep="/") )
source( paste(mchanlcdir, "scriptR/fista.lib.R", sep="/") )
options(digits=10)
library(plotrix)

fista.mat2.weight.cv <- function(infile1, infile2, freq.file, lambda.file, nfold, outdir)
{
    printf(" --- fista.mat2.weight.cv ---- \n")
    printf(" infile1     = %s\n", infile1)
    printf(" infile2     = %s\n", infile2)
    printf(" freq.file   = %s\n", freq.file)
    printf(" lambda.file = %s\n", lambda.file)
    printf(" nfold       = %d\n", nfold)
    printf(" outdir      = %s\n", outdir)
    printf(" --- fista.mat2.weight.cv ---- \n")

    dir.create(outdir)

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
    ## add random index
    ##
    ri.vec1 = GetRandIndex(nrow1)
    ri.vec2 = GetRandIndex(nrow2)
    ## data1.ri.df = cbind(data1.df, ri.vec1)
    ## data2.ri.df = cbind(data2.df, ri.vec2)


    ##
    ## read lambda file
    ##
    lambda.df <- read.table(lambda.file, comment.char = "#")
    lambda.lo    = lambda.df[1, 1]
    lambda.up    = lambda.df[1, 2]
    nlambda      = lambda.df[1, 3]
    lambda.vec   = GetPointLog(nlambda, lambda.lo, lambda.up)

    print(lambda.vec)

    rms.mean.vec = numeric(0)
    rms.sd.vec = numeric(0)
   
    for(ilambda in 1:nlambda){
        lambda = lambda.vec[ilambda]
        printf("lambda = %e\n", lambda)

        outdir.lambda = sprintf("%s/lambda_%2.2d", outdir, ilambda)
        dir.create(outdir.lambda)

        ##
        ## Get Nfold Index
        ##
        nterm = nfold
        nelem1 = as.integer(nrow1 / nterm)
        nelem2 = as.integer(nrow2 / nterm)

        rms.vec = numeric(0)
        for(iterm in 1:nterm){
            outdir.term = sprintf("%s/cv_%2.2d", outdir.lambda, iterm)
            dir.create(outdir.term)
            
            index1.st = (iterm - 1) * nelem1 + 1
            index1.ed = iterm * nelem1
            if(iterm == nterm){
                index1.ed = nrow1
            }
            data1.tr.df = data1.df[ri.vec1[-(index1.st:index1.ed)],]
            data1.te.df = data1.df[ri.vec1[index1.st:index1.ed],]

            index2.st = (iterm - 1) * nelem2 + 1
            index2.ed = iterm * nelem2
            if(iterm == nterm){
                index2.ed = nrow2
            }
            data2.tr.df = data2.df[ri.vec2[-(index2.st:index2.ed)],]
            data2.te.df = data2.df[ri.vec2[index2.st:index2.ed],]

            ##
            ## sort
            ##
            data1.tr.sort.df = data1.tr.df[order(data1.tr.df[,1]),]
            data1.te.sort.df = data1.te.df[order(data1.te.df[,1]),]
            data2.tr.sort.df = data2.tr.df[order(data2.tr.df[,1]),]
            data2.te.sort.df = data2.te.df[order(data2.te.df[,1]),]

            ##
            ## write df
            ##
            data1.tr.file = sprintf("%s/data1_tr.txt", outdir.term)
            data1.te.file = sprintf("%s/data1_te.txt", outdir.term)
            data2.tr.file = sprintf("%s/data2_tr.txt", outdir.term)
            data2.te.file = sprintf("%s/data2_te.txt", outdir.term)
            
            write.table(data1.tr.sort.df, data1.tr.file, row.names = FALSE, col.names = FALSE)
            write.table(data1.te.sort.df, data1.te.file, row.names = FALSE, col.names = FALSE)
            write.table(data2.tr.sort.df, data2.tr.file, row.names = FALSE, col.names = FALSE)
            write.table(data2.te.sort.df, data2.te.file, row.names = FALSE, col.names = FALSE)

            ##
            ## for training data
            ##
            
            ##
            ## make fourier matrix file
            ##
            cmd.bin = "/home/morii/work/github/moriiism/mchanlc/mkmat/mkmat2"
            cmd = sprintf("%s  %s  %s  %s  %s  %s",
                cmd.bin, data1.tr.file, data2.tr.file, freq.file, outdir.term, "mat")
            system(cmd)

            ##
            ## read fourier matrix file
            ##
            mat.w.file = sprintf("%s/mat_12_w.dat", outdir.term)
            mat.w.df <- read.table(mat.w.file, comment.char = "#")
            A.w.mat = data.matrix(mat.w.df)
            printf("nrow(A.w.mat) = %d\n", nrow(A.w.mat))
            printf("ncol(A.w.mat) = %d\n", ncol(A.w.mat))

            mat.file = sprintf("%s/mat_12.dat", outdir.term)
            mat.df <- read.table(mat.file, comment.char = "#")
            A.mat = data.matrix(mat.df)


            ##
            ## for test data
            ##

            ##
            ## make fourier matrix file
            ##
            cmd.bin = "/home/morii/work/github/moriiism/mchanlc/mkmat/mkmat2"
            cmd = sprintf("%s  %s  %s  %s  %s  %s",
                cmd.bin, data1.te.file, data2.te.file, freq.file, outdir.term, "mat_te")
            system(cmd)

            ##
            ## read fourier matrix file
            ##
            mat.te.w.file = sprintf("%s/mat_te_12_w.dat", outdir.term)
            mat.te.w.df <- read.table(mat.te.w.file, comment.char = "#")
            A.w.mat.te = data.matrix(mat.te.w.df)
            printf("nrow(A.w.mat.te) = %d\n", nrow(A.w.mat.te))
            printf("ncol(A.w.mat.te) = %d\n", ncol(A.w.mat.te))

            mat.te.file = sprintf("%s/mat_te_12.dat", outdir.term)
            mat.te.df <- read.table(mat.te.file, comment.char = "#")
            A.mat.te = data.matrix(mat.te.df)
            
            rms = fista.mat2.weight(data1.tr.sort.df, data2.tr.sort.df, freq.file,
                A.mat, A.w.mat, lambda, outdir.term,
                data1.te.sort.df, data2.te.sort.df, A.mat.te)
            rms.vec = c(rms.vec, rms)
        }
        rms.mean.vec = c(rms.mean.vec, mean(rms.vec))
        rms.sd.vec = c(rms.sd.vec, sd(rms.vec))

        printf("-----\n")
        print(rms.mean.vec)
        print(rms.sd.vec)
        printf("-----\n")
    }

    cv.df = cbind(lambda.vec, rms.mean.vec, rms.sd.vec)
    print(lambda.vec)
    print(rms.mean.vec)
    print(rms.sd.vec)
    
    cv.file = sprintf("%s/cv.dat", outdir)
    write(t(cv.df), cv.file, ncolumns = 3)
}


fista.mat2.weight <- function(data1.df, data2.df, freq.file, A.mat, A.w.mat, lambda, outdir, data1.te.df, data2.te.df, A.mat.te)
{
    ##
    ## read frequency setup file
    ##
    freq.df <- read.table(freq.file, comment.char = "#")
    freq.lo = freq.df[1, 1]
    freq.up = freq.df[1, 2]
    nfreq   = freq.df[1, 3]
    delta.freq = (freq.up - freq.lo) / nfreq
    ### ncol = nfreq * 2 * 2
    
    nrow1 = nrow(data1.df)
    nrow2 = nrow(data2.df)
    ncol  = ncol(A.mat)

    data.df = rbind(data1.df, data2.df)
    
    h.vec = as.vector(c(data1.df[,2] * sqrt(nrow2), data2.df[,2] * sqrt(nrow1) ))

    hist.lo = -0.005
    hist.up =  0.005
    hist.nbin = 200
    hist.delta = (hist.up - hist.lo) / hist.nbin
    
    tolerance = 5.e-8
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
        ik = FindIk(y, L.pre, eta, A.w.mat, h.vec, lambda)
        ## printf("ik = %d\n", ik)
        L = eta**ik * L.pre
        x = ProxMap(y, L, A.w.mat, h.vec, lambda)
        t.new = (1. + sqrt(1. + 4 * t**2))/2.
        y.new = x + (t - 1.)/ t.new * (x - x.pre)
        x.pre = x
        y = y.new
        t = t.new
        L.pre = L

        cost = FGFunc(y.new, A.w.mat, h.vec, lambda)
        ##printf("ik = %d, L = %e, t = %e, FGFunc = %e\n",
        ##       ik, L, t, FGFunc(y.new, A.w.mat, h.vec, lambda))

        ##printf("diff = %e\n", (cost.pre - cost) / cost)
        
        if(k > 1 && (cost.pre - cost) / cost < tolerance){
            printf("k = %d, cost = %e\n", k, cost)
            break
        }
        cost.pre = cost
    }

###    n.nonzero[ilambda] = 0
###    for(i in 1:length(x)){
###        if(abs(x[i]) > lambda){
###            n.nonzero[ilambda] = n.nonzero[ilambda] + 1
###        }
###    }
###    printf("n.nonzero = %d\n", n.nonzero[ilambda])

    ##
    ## delta.time
    ##
    N = ncol/4
    delta.time12.vec = numeric(0)
    weight.delta.time12.vec = numeric(0)

    delta.time12.vec.tmp = numeric(0)
    weight.delta.time12.vec.tmp = numeric(0)
    period.vec.tmp = numeric(0)
    
    delta.time1.vec = numeric(0)
    delta.time2.vec = numeric(0)
    isel.vec = numeric(0)

    x.common = rep(0.0, ncol)
    
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

            ## select common between X and Opt
            if(0.0 < delta.time12 && delta.time12 < 0.01){
                x.common[i] = x[i]
                x.common[i + N] = x[i + N]
                x.common[i + 2 * N] = x[i + 2 * N]
                x.common[i + 3 * N] = x[i + 3 * N]
            }
            
            ##printf("delta.time1 = %e\n", delta.time1)
            ##printf("delta.time2 = %e\n", delta.time2)                
            ##printf("delta.time12 = %e\n", delta.time12)
            ##printf("T = %e\n", 1./freq)
            ##printf("weight.delta.time12 = %e\n", weight.delta.time12)

            delta.time12.vec.tmp = c(delta.time12.vec.tmp, delta.time1 - delta.time2)
            weight.delta.time12.vec.tmp = c(weight.delta.time12.vec.tmp, weight.delta.time12)
            period.vec.tmp = c(period.vec.tmp, 1./freq)
            
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

    ##
    ## lag v.s. norm
    ##
    lag.norm = cbind(delta.time12.vec.tmp, weight.delta.time12.vec.tmp)
    outfile = sprintf("%s/lag_norm.dat", outdir)
    write(t(lag.norm), file=outfile, ncolumns = 2)

    ##
    ## lag v.s. norm v.s. period
    ##
    lag.norm = cbind(delta.time12.vec.tmp, weight.delta.time12.vec.tmp, period.vec.tmp)
    outfile = sprintf("%s/lag_norm_period.dat", outdir)
    write(t(lag.norm), file=outfile, ncolumns = 3)

    h.rec.vec = A.mat %*% x

    ##
    ## h.rec.vec = h.rec.vec * sd + mean
    ##
    
    lc.rec = cbind(data.df[, 1], h.rec.vec)
    outfile = sprintf("%s/rec.dat", outdir)
    write(t(lc.rec), file=outfile, ncolumns = 2)

    outfile = sprintf("%s/x.dat", outdir)
    x.mat = cbind(c(1:(length(x)/2),1:(length(x)/2)), x)
    write(t(x.mat), file=outfile, ncolumns = 2)

    ##
    ## residual for test data
    ##
    data.te.df = rbind(data1.te.df, data2.te.df)
    h.rec.te.vec = A.mat.te %*% x
    rms = sqrt( sum( (data.te.df[, 2] - h.rec.te.vec) * (data.te.df[, 2] - h.rec.te.vec) ) )
    outfile = sprintf("%s/rms.dat", outdir)
    write(rms, file=outfile)
    
    ##
    ## common
    ##
    h.rec.common.vec = A.mat %*% x.common

    lc.common.rec = cbind(data.df[, 1], h.rec.common.vec)
    outfile = sprintf("%s/rec_common.dat", outdir)
    write(t(lc.common.rec), file=outfile, ncolumns = 2)

    outfile = sprintf("%s/x_common.dat", outdir)
    x.common.mat = cbind(c(1:(length(x.common)/2),1:(length(x.common)/2)), x.common)
    write(t(x.common.mat), file=outfile, ncolumns = 2)

    power.common = 0.0
    for(i in 1:N){
        power.common = power.common + x.common[i]**2 + x.common[i + N]**2 + x.common[i + 2 * N]**2 + x.common[i + 3 * N]**2
    }
    printf("power.common = %e\n", power.common)

    power = 0.0
    for(i in 1:N){
        power = power + x[i]**2 + x[i + N]**2 + x[i + 2 * N]**2 + x[i + 3 * N]**2
    }
    printf("power = %e\n", power)

    
    ##
    ## uncommon 
    ##
    
    lc.uncommon.rec = cbind(data.df[, 1], h.rec.vec - h.rec.common.vec)
    outfile = sprintf("%s/rec_uncommon.dat", outdir)
    write(t(lc.uncommon.rec), file=outfile, ncolumns = 2)

    return(rms)
}
