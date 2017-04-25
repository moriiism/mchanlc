###
### grplasso.R
###
### 2017.03.30 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
options(digits=10)
library(grplasso)

GrpLasso <- function(infile, freq.file, mat.file, nlambda)
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

    ##
    ## calc lambda.array
    ##
    B.mat = cbind(1, A.mat)
    h.vec = as.vector(data.df[,2])
    index.vec = c(NA, 1:(ncol/2), 1:(ncol/2))

    lambda.array <- lambdamax(B.mat, y = h.vec, index = index.vec,
                              penscale = sqrt, model = LinReg()) * 0.5^(0:(nlambda-2))
    lambda.array = c(lambda.array, 0.0)
    
###    fit = grplasso(x = B.mat, y = h.vec, index = index.vec, model = LinReg(),
###        lambda = lambda.array, penscale = sqrt, center = TRUE,
###        standardize = TRUE,
###        control = grpl.control(update.hess = "lambda", trace = 0))

    for(ilambda in 1:nlambda){

        fit = grplasso(x = B.mat, y = h.vec, index = index.vec, model = LinReg(),
            lambda = lambda.array[ilambda], center = FALSE,
            standardize = FALSE, control = grpl.control(update.hess = "lambda", update.every = 1))
        print(fit)
        write(coef(fit), file="fit.dat")

        res.vec = coef(fit)[,1]

        h.rec.vec = A.mat %*% res.vec[2:(ncol+1)]
        lc.rec = cbind(data.df[, 1], h.rec.vec)
        outfile = sprintf("rec_grplasso_%2.2d.dat", ilambda)
        write(t(lc.rec), file=outfile, ncolumns = 2)

        outfile = sprintf("x_grplasso_%2.2d.dat", ilambda)
        res.mat = cbind(c(1:(ncol/2),1:(ncol/2)), res.vec[2:(ncol+1)])
        write(t(res.mat), file=outfile, ncolumns = 2)
    }

    for(ilambda in 1:nlambda){
        lambda = lambda.array[ilambda]
        printf("lambda = %e\n", lambda)
    }
}
