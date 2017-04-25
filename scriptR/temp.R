mitooldir = "/home/morii/work/github/moriiism/mitool"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
options(digits=10)

### ibin = 1, ..., nbin

GetIbin <- function(xval, xlo, xup, nbinx)
{
    delta.xval = (xup - xlo) / nbinx
    ibin = ceiling((xval - xlo) / delta.xval)
    return(ibin)
}

PutHist <- function(xval, weight, hist.xval, hist.oval)
{
    ibin = GetIbin(xval, xlo, xup, nbin)
    hist.oval[ibin] = hist.oval[ibin] + weight
}

PlotHist <- function(x, lambda)
{
    


    
}



