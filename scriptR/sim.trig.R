###
### sim.trig.R
###
### Simulate light curve of linear combination of
### trigonometric functions (sine function),
### using frequency and phase information, on GTI duration.
###
### 2017.03.31 M.Morii
###
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"

###

source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
options(digits=10)

TrigFunc <- function(time, trig.df)
{
    nrow = nrow(trig.df)
    ans = 0.0
    for(irow in 1:nrow){
        freq  = as.numeric(trig.df[irow, 1])
        phase = as.numeric(trig.df[irow, 2])
        norm  = as.numeric(trig.df[irow, 3])
        lag   = as.numeric(trig.df[irow, 4])
        ans = ans + norm * sin(2 * pi * (freq * (time - lag)  - phase))
    }
    return(ans)
}

GtiIn <- function(time, gti.df)
{
    nrow = nrow(gti.df)
    ans = 0
    for(irow in 1:nrow){
        if(gti.df[irow, 1] <= time && time <= gti.df[irow, 2]){
            ans = 1
            break
        }
    }
    return(ans)
}

GtiItime <- function(time.lo, delta.time, gti.df)
{
    gti.itime.df = gti.df
    nrow = nrow(gti.df)
    ## start time
    for(irow in 1:nrow){
        gti.itime.df[irow, 1] = as.integer(ceiling((gti.df[irow, 1] - time.lo) / delta.time ))
    }
    ## end time
    for(irow in 1:nrow){
        gti.itime.df[irow, 2] = as.integer(floor((gti.df[irow, 2] - time.lo) / delta.time ))
    }
    return(gti.itime.df)
}

