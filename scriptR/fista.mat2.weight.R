###
### fista.mat2.R
###
### 2017.03.30 M.Morii
### 2017.04.19 M.Morii
###  for 2 light curve
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlcdir = "/home/morii/work/github/moriiism/mchanlc"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mitooldir, "script/mirlib/binning.R", sep="/") )
source( paste(mchanlcdir, "scriptR/fista.lib.R", sep="/") )
options(digits=10)
library(plotrix)

fista.mat2.weight <- function(infile1, infile2, freq.file, lambda, outdir)
{
    printf(" --- fista.mat2.weight ---- \n")
    printf(" infile1    = %s\n", infile1)
    printf(" infile2    = %s\n", infile2)
    printf(" freq.file  = %s\n", freq.file)
    printf(" lambda     = %e\n", lambda)
    printf(" outdir     = %s\n", outdir)
    printf(" --- fista.mat2.weight ---- \n")

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
    ## make fourier matrix file
    ##
    cmd.bin = "/home/morii/work/github/moriiism/mchanlc/mkmat/mkmat2"
    cmd = sprintf("%s  %s  %s  %s  %s  %s",
        cmd.bin, infile1, infile2, freq.file, outdir, "mat")
    system(cmd)

    ##
    ## read fourier matrix file
    ##
    mat.w.file = sprintf("%s/mat_12_w.dat", outdir)
    mat.w.df <- read.table(mat.w.file, comment.char = "#")
    A.w.mat = data.matrix(mat.w.df)
    printf("nrow(A.w.mat) = %d\n", nrow(A.w.mat))
    printf("ncol(A.w.mat) = %d\n", ncol(A.w.mat))

    mat.file = sprintf("%s/mat_12.dat", outdir)
    mat.df <- read.table(mat.file, comment.char = "#")
    A.mat = data.matrix(mat.df)

#######

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

    ##
    ## delta.time
    ##
    N = ncol/4
    delta.time12.vec = numeric(0)
    norm.delta.time12.vec = numeric(0)
    norm.opt.delta.time12.vec = numeric(0)
    norm.xray.delta.time12.vec = numeric(0)
    period.vec = numeric(0)

    x.lagplus  = rep(0.0, ncol)
    x.lagminus = rep(0.0, ncol)
    
    for(i in 1:N){
        if(abs(x[i]) > 0.0  && abs(x[i + N]) > 0.0 && abs(x[i + 2 * N]) > 0.0 && abs(x[i + 3 * N]) > 0.0 ) {
            freq = freq.lo + delta.freq * (i - 0.5)
            omega = 2 * pi * freq
            delta.time1 = GetAlphaOfSinCos(x[i + N], x[i]) / omega
            delta.time2 = GetAlphaOfSinCos(x[i + 3 * N], x[i + 2 * N]) / omega
            delta.time1 = GetDeltaTime(delta.time1, 1./freq)
            delta.time2 = GetDeltaTime(delta.time2, 1./freq)
            delta.time12 = delta.time1 - delta.time2
            
            norm.delta.time12      = sqrt( x[i]**2 + x[i + N]**2 + x[i + 2 * N]**2 + x[i + 3 * N]**2 )
            norm.opt.delta.time12  = sqrt( x[i]**2 + x[i + N]**2)
            norm.xray.delta.time12 = sqrt( x[i + 2 * N]**2 + x[i + 3 * N]**2)

            ## select lagplus between X and Opt
            if(0.0 <= delta.time12 ){
                x.lagplus[i] = x[i]
                x.lagplus[i + N] = x[i + N]
                x.lagplus[i + 2 * N] = x[i + 2 * N]
                x.lagplus[i + 3 * N] = x[i + 3 * N]
            }
            if(0.0 > delta.time12 ){
                x.lagminus[i] = x[i]
                x.lagminus[i + N] = x[i + N]
                x.lagminus[i + 2 * N] = x[i + 2 * N]
                x.lagminus[i + 3 * N] = x[i + 3 * N]
            }


            
            delta.time12.vec = c(delta.time12.vec, delta.time12)
            norm.delta.time12.vec = c(norm.delta.time12.vec, norm.delta.time12)
            norm.opt.delta.time12.vec = c(norm.opt.delta.time12.vec, norm.opt.delta.time12)
            norm.xray.delta.time12.vec = c(norm.xray.delta.time12.vec, norm.xray.delta.time12)
            period.vec = c(period.vec, 1./freq)
        }
    }

    ##
    ## lag, norm, norm.opt, norm.xray, period
    ##
    sum.df = cbind(delta.time12.vec, norm.delta.time12.vec, norm.opt.delta.time12.vec, norm.xray.delta.time12.vec, period.vec)
    outfile = sprintf("%s/sum.dat", outdir)
    write(t(sum.df), file=outfile, ncolumns = 5)

    h.rec.vec = A.mat %*% x

    lc.rec = cbind(data.df[, 1], h.rec.vec)
    outfile = sprintf("%s/rec.dat", outdir)
    write(t(lc.rec), file=outfile, ncolumns = 2)

    outfile = sprintf("%s/x.dat", outdir)
    x.mat = cbind(c(1:(length(x)/2),1:(length(x)/2)), x)
    write(t(x.mat), file=outfile, ncolumns = 2)


    ##
    ## lagplus
    ##
    h.rec.lagplus.vec = A.mat %*% x.lagplus

    lc.lagplus.rec = cbind(data.df[, 1], h.rec.lagplus.vec)
    outfile = sprintf("%s/rec_lagplus.dat", outdir)
    write(t(lc.lagplus.rec), file=outfile, ncolumns = 2)

    outfile = sprintf("%s/x_lagplus.dat", outdir)
    x.lagplus.mat = cbind(c(1:(length(x.lagplus)/2),1:(length(x.lagplus)/2)), x.lagplus)
    write(t(x.lagplus.mat), file=outfile, ncolumns = 2)


    ##
    ## lagminus
    ##
    h.rec.lagminus.vec = A.mat %*% x.lagminus

    lc.lagminus.rec = cbind(data.df[, 1], h.rec.lagminus.vec)
    outfile = sprintf("%s/rec_lagminus.dat", outdir)
    write(t(lc.lagminus.rec), file=outfile, ncolumns = 2)

    outfile = sprintf("%s/x_lagminus.dat", outdir)
    x.lagminus.mat = cbind(c(1:(length(x.lagminus)/2),1:(length(x.lagminus)/2)), x.lagminus)
    write(t(x.lagminus.mat), file=outfile, ncolumns = 2)


    ##
    ## lagplus + lagminus
    ##
    h.rec.lagboth.vec =  h.rec.lagminus.vec + h.rec.lagplus.vec
    lc.lagboth.rec = cbind(data.df[, 1], h.rec.lagboth.vec)
    outfile = sprintf("%s/rec_lagboth.dat", outdir)
    write(t(lc.lagboth.rec), file=outfile, ncolumns = 2)

    
    ##
    ## power
    ##
    power.opt.lagplus = 0.0
    power.xray.lagplus = 0.0
    power.opt.lagminus = 0.0
    power.xray.lagminus = 0.0    
    for(i in 1:N){
        power.opt.lagplus   = power.opt.lagplus + x.lagplus[i]**2 + x.lagplus[i + N]**2
        power.xray.lagplus  = power.xray.lagplus + x.lagplus[i + 2 * N]**2 + x.lagplus[i + 3 * N]**2
        power.opt.lagminus  = power.opt.lagminus + x.lagminus[i]**2 + x.lagminus[i + N]**2
        power.xray.lagminus = power.xray.lagminus + x.lagminus[i + 2 * N]**2 + x.lagminus[i + 3 * N]**2        
    }

    printf("power.opt.lagplus  = %e\n", power.opt.lagplus)
    printf("power.opt.lagminus = %e\n", power.opt.lagminus)
    printf("power.xray.lagplus  = %e\n", power.xray.lagplus)
    printf("power.xray.lagminus = %e\n", power.xray.lagminus)        
}
    
