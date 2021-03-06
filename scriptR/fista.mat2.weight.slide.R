###
### fista.mat2.weight.slide.R
###
### 2017.05.23 M.Morii
###  sliding window
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
mchanlcdir = "/home/morii/work/github/moriiism/mchanlc"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(mitooldir, "script/mirlib/binning.R", sep="/") )
source( paste(mchanlcdir, "scriptR/fista.lib.R", sep="/") )
options(digits=10)
library(plotrix)


GetDataInWin <- function(data.df, win.lo, win.up)
{
    data.sel.df = data.df[win.lo <= data.df$V1 & data.df$V1 <= win.up,]
    return(data.sel.df)
}


fista.mat2.weight.slide <- function(infile1, infile2, freq.file, lambda, slide.file, outdir)
{
    printf(" --- fista.mat2.weight ---- \n")
    printf(" infile1    = %s\n", infile1)
    printf(" infile2    = %s\n", infile2)
    printf(" freq.file  = %s\n", freq.file)
    printf(" lambda     = %e\n", lambda)
    printf(" slide.file = %s\n", slide.file)
    printf(" outdir     = %s\n", outdir)
    printf(" --- fista.mat2.weight ---- \n")

    dir.create(outdir)

    ##
    ## read slide file
    ##
    slide.df <- read.table(slide.file, comment.char = "#")
    slide.lo     = slide.df[1, 1]
    slide.up     = slide.df[1, 2]
    slide.npoint = slide.df[1, 3]
    slide.width  = slide.df[1, 4]

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
    ## loop
    ##
    slide.vec  = GetPointLin(slide.npoint, slide.lo, slide.up)

    power.opt.lagplus.vec    = numeric(0)
    power.opt.lagminus.vec   = numeric(0)
    power.xray.lagplus.vec   = numeric(0)
    power.xray.lagminus.vec  = numeric(0)
    ratio.lag.power.opt.vec  = numeric(0)
    ratio.lag.power.xray.vec = numeric(0)
    for(islide in 1:slide.npoint){
        win.lo = slide.vec[islide]
        win.up = win.lo + slide.width
        outdir.slide = sprintf("%s/slide_%2.2d", outdir, islide)
        dir.create(outdir.slide)

        slide.center = win.lo + slide.width / 2.0

        ##
        ## select data
        ##
        data1.sel.df = GetDataInWin(data1.df, win.lo, win.up)
        data2.sel.df = GetDataInWin(data2.df, win.lo, win.up)
        outfile1.sel = sprintf("%s/data1.dat", outdir.slide)
        write(t(data1.sel.df), file=outfile1.sel, ncolumns = 2)
        outfile2.sel = sprintf("%s/data2.dat", outdir.slide)
        write(t(data2.sel.df), file=outfile2.sel, ncolumns = 2)

        power.ret = FistaMat2Weight(outfile1.sel, outfile2.sel, freq.file, lambda, slide.width, slide.center, outdir.slide)

        power.opt.lagplus   = power.ret[1]
        power.opt.lagminus  = power.ret[2]
        power.xray.lagplus  = power.ret[3]
        power.xray.lagminus = power.ret[4]
        power.opt.lagplus.vec   = c(power.opt.lagplus.vec, power.opt.lagplus)
        power.opt.lagminus.vec  = c(power.opt.lagminus.vec, power.opt.lagminus)
        power.xray.lagplus.vec  = c(power.xray.lagplus.vec, power.xray.lagplus)
        power.xray.lagminus.vec = c(power.xray.lagminus.vec, power.xray.lagminus)

        ratio.lag.power.opt  = power.opt.lagplus / (power.opt.lagplus + power.opt.lagminus)
        ratio.lag.power.xray = power.xray.lagplus / (power.xray.lagplus + power.xray.lagminus)
        ratio.lag.power.opt.vec = c(ratio.lag.power.opt.vec, ratio.lag.power.opt)
        ratio.lag.power.xray.vec = c(ratio.lag.power.xray.vec, ratio.lag.power.xray)
        
    }

    slide.center.vec = slide.vec + slide.width / 2.0
    power.df = cbind(slide.center.vec, power.opt.lagplus.vec, power.opt.lagminus.vec, power.xray.lagplus.vec, power.xray.lagminus.vec)
    outfile = sprintf("%s/pow.dat", outdir)
    write(t(power.df), file=outfile)

    WritePowQdp(data1.df, data2.df, slide.center.vec,
                power.opt.lagplus.vec, power.xray.lagplus.vec, power.opt.lagminus.vec, power.xray.lagminus.vec, outdir)
    WriteRatioLagQdp(data1.df, data2.df, slide.center.vec,
                     ratio.lag.power.opt.vec, ratio.lag.power.xray.vec, outdir)
}

WritePowQdp <- function(data1.df, data2.df, slide.center.vec,
                        power.opt.lagplus.vec, power.xray.lagplus.vec, power.opt.lagminus.vec, power.xray.lagminus.vec, outdir)
{
    nrow1 = nrow(data1.df)
    nrow2 = nrow(data2.df)

    
    outfile = sprintf("%s/pow.qdp", outdir)
    sink(outfile)
    printf("skip sing\n")
    printf("\n")
    printf("! optical\n")
    for(irow in 1:nrow1){
        printf("%e %e\n", data1.df[irow, 1], data1.df[irow, 2])
    }
    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! X-ray\n")
    for(irow in 1:nrow2){
        printf("%e %e\n", data2.df[irow, 1], data2.df[irow, 2])
    }

    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! power opt.lagplus\n")
    for(inum in 1:length(slide.center.vec)){
        printf("%e %e\n", slide.center.vec[inum], power.opt.lagplus.vec[inum])
    }
    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! power xray.lagplus\n")
    for(inum in 1:length(slide.center.vec)){
        printf("%e %e\n", slide.center.vec[inum], power.xray.lagplus.vec[inum])
    }

    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! power opt.lagminus\n")
    for(inum in 1:length(slide.center.vec)){
        printf("%e %e\n", slide.center.vec[inum], power.opt.lagminus.vec[inum])
    }
    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! power xray.lagminus\n")
    for(inum in 1:length(slide.center.vec)){
        printf("%e %e\n", slide.center.vec[inum], power.xray.lagminus.vec[inum])
    }
    printf("\n")

    printf("la file\n")
    printf("time off\n")
    printf("lw 5\n")
    printf("csize 1.2\n")
    printf("la rot\n")
    printf("\n")
    printf("win 1\n")
    printf("yplot 1 2\n")
    printf("loc  0.1 0.44999999 1 0.94999999\n")    
    printf("la y Light Curve\n")
    printf("la pos y 5.0\n")
    printf("lab nx off\n")
    printf("\n")
    printf("win 2\n")
    printf("yplot 3 4 5 6\n")
    printf("loc  0.1 4.99999821E-2 1 0.54999995\n")
    printf("la y Power\n")
    printf("la pos y 5.0\n")

    printf("col 1 on 3\n")
    printf("col 2 on 4\n")
    printf("lst 1 on 3\n")
    printf("lst 1 on 4\n")    

    printf("col 1 on 5\n")
    printf("col 2 on 6\n")
    printf("lst 2 on 5\n")
    printf("lst 2 on 6\n")    
  
    printf("la x time (MJD - 57388)\n")
    printf("win all\n")
    printf("\n")
    printf("\n")
    printf("\n")    
    sink()
}    


WriteRatioLagQdp <- function(data1.df, data2.df, slide.center.vec,
                             ratio.lag.power.opt.vec, ratio.lag.power.xray.vec, outdir)

{
    nrow1 = nrow(data1.df)
    nrow2 = nrow(data2.df)
    
    outfile = sprintf("%s/ratio_lag_pm.qdp", outdir)
    sink(outfile)
    printf("skip sing\n")
    printf("\n")
    printf("! optical\n")
    for(irow in 1:nrow1){
        printf("%e %e\n", data1.df[irow, 1], data1.df[irow, 2])
    }
    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! X-ray\n")
    for(irow in 1:nrow2){
        printf("%e %e\n", data2.df[irow, 1], data2.df[irow, 2])
    }

    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! ratio.lag.power.opt\n")
    for(inum in 1:length(slide.center.vec)){
        printf("%e %e\n", slide.center.vec[inum], ratio.lag.power.opt.vec[inum])
    }
    printf("\n")
    printf("no\n")
    printf("\n")
    printf("! ratio.lag.power.xray\n")
    for(inum in 1:length(slide.center.vec)){
        printf("%e %e\n", slide.center.vec[inum], ratio.lag.power.xray.vec[inum])
    }
    printf("\n")

    printf("la file\n")
    printf("time off\n")
    printf("lw 5\n")
    printf("csize 1.2\n")
    printf("la rot\n")
    printf("\n")
    printf("win 1\n")
    printf("yplot 1 2\n")
    printf("loc  0.1 0.44999999 1 0.94999999\n")    
    printf("la y Light Curve\n")
    printf("la pos y 5.0\n")
    printf("lab nx off\n")
    printf("\n")
    printf("win 2\n")
    printf("yplot 3 4 \n")
    printf("loc  0.1 4.99999821E-2 1 0.54999995\n")
    printf("la y P(+) / [P(+) + P(-)] \n")
    printf("la pos y 5.0\n")
    printf("r y 0 1 \n")

    printf("col 1 on 3\n")
    printf("col 2 on 4\n")
  
    printf("la x time (MJD - 57388)\n")
    printf("win all\n")
    printf("\n")
    printf("\n")
    printf("\n")    
    sink()
}    







FistaMat2Weight <-function(infile1, infile2, freq.file, lambda, slide.width, slide.center, outdir)
{
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
### ncol = nfreq * 2 * 2
        
    nrow1 = nrow(data1.df)
    nrow2 = nrow(data2.df)
    ncol  = ncol(A.mat)

    data.df = rbind(data1.df, data2.df)
        
    h.vec = as.vector(c(data1.df[,2] * sqrt(nrow2), data2.df[,2] * sqrt(nrow1) ))
        
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

            period = 1./freq
            if(period < slide.width){
                
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

    sum.2.df = cbind(slide.center, delta.time12.vec, norm.delta.time12.vec, norm.opt.delta.time12.vec, norm.xray.delta.time12.vec, period.vec)
    outfile = sprintf("%s/sum2.dat", outdir)
    write(t(sum.2.df), file=outfile, ncolumns = 6)    

    
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

    outfile = sprintf("%s/power.dat", outdir)
    sink(outfile)
    printf("%e ! power.opt.lagplus\n", power.opt.lagplus)
    printf("%e ! power.opt.lagminus\n", power.opt.lagminus)
    printf("%e ! power.xray.lagplus\n", power.xray.lagplus)
    printf("%e ! power.xray.lagminus\n", power.xray.lagminus)
    sink()

    list = c(power.opt.lagplus, power.opt.lagminus, power.xray.lagplus, power.xray.lagminus)
    return(list)
    
}









    
