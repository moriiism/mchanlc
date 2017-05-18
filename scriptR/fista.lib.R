###
### fista.lib.R
###
### 2017.05.18 M.Morii
###   library for fista
###

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
