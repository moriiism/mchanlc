###
### fista.mat.R
###
### 2017.03.30 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
options(digits=10)
library(grplasso)

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

SoftThres2 <- function(a, b, L, lambda)
{
    l = sqrt(a**2 + b**2) * L
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

fista.mat <- function(infile, freq.file, mat.file, nlambda)
{
    printf(" --- fista.mat ---- \n")
    printf(" infile    = %s\n", infile)
    printf(" freq.file = %s\n", freq.file)
    printf(" mat.file  = %s\n", mat.file)
    printf(" nlambda   = %d\n", nlambda)
    printf(" --- fista.mat ---- \n")

    ##
    ## read infile
    ##
    data.df <- read.table(infile, comment.char = "#")
    printf("ncol(data.df) = %d\n", ncol(data.df))
    printf("nrow(data.df) = %d\n", nrow(data.df))
    nrow = nrow(data.df)

    ##
    ## read frequency setup file
    ##
    freq.df <- read.table(freq.file, comment.char = "#")
    freq.lo = freq.df[1, 1]
    freq.up = freq.df[1, 2]
    nfreq   = freq.df[1, 3]

    delta.freq = (freq.up - freq.lo) / nfreq
    ncol = nfreq * 2

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
    mean = mean(h.vec)
    sd   = sd(h.vec)
    h.vec = h.vec - mean
    h.vec = h.vec / sd
    
    ##
    ## calc lambda.array
    ##
    index.vec = c(NA, 1:(ncol/2), 1:(ncol/2))

    lambda.array <- lambdamax(B.mat, y = h.vec, index = index.vec,
                              penscale = sqrt, model = LinReg()) * 0.5^(0:(nlambda-2))
    lambda.array = c(lambda.array, 0.0)

    #### ad hoc
    lambda.array = 1e3 * 0.5^(0:(nlambda-2))
    lambda.array = c(lambda.array, 0.0)
    ####

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
        
        h.rec.vec = A.mat %*% x

        ##
        h.rec.vec = h.rec.vec * sd + mean
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



