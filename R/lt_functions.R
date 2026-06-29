# For classic abridged data (i.e. for ages <1, 1-4, 5-9, ..., 90+)
lt_abridged <- function (nMx, x = c(0, 1, cumsum(rep(5, length(nMx) - 2)))) {
  N <- length(nMx)
  w <- c(diff(x), 5)
  # You may what to change the ax for <1 and 1.4 
  # also these are not sex specific
  ax <- c(0.07 + 1.7 * nMx[1], 0.4, rep(0.5, N - 3), 1/nMx[N]) 
  nax <- w * ax
  nqx <- (w * nMx)/(1 + w * (1 - ax) * nMx)
  nqx[N] <- 1
  npx <- 1 - nqx
  lx <- c(1, cumprod(npx[-N]))
  ndx <- -diff(c(lx, 0))
  nLx <- lx[-1] * w[-N] + ndx[-N] * nax[-N]
  nLx[N] <- ax[N] * lx[N]
  Tx <- rev(cumsum(rev(nLx)))
  ex <- Tx/lx
  
  data.frame(
    x=x,nMx=nMx,nax=nax,
    nqx=nqx,npx=npx,ndx=ndx,
    lx=lx,nLx=nLx,Tx=Tx,ex=ex
    )
  }

# For single year data (i.e. for ages 0, 1, 2, ..., 100+)
lt_single <- function (Mx, x = 0:(length(Mx)-1)) {
  N <- length(Mx)
  ax <- c(0.14, rep(0.5, N - 2), 1/Mx[N]) # You may what to change the ax for <1
  qx <- Mx/(1 + (1 - ax) * Mx)
  qx[N] <- 1
  px <- 1 - qx
  lx <- c(1, cumprod(px[-N]))
  dx <- -diff(c(lx, 0))
  Lx <- lx[-1] + dx[-N] * ax[-N]
  Lx[N] <- ax[N] * lx[N]
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/lx
  
  data.frame(
    x=x,Mx=Mx,ax=ax,
    qx=qx,px=px,dx=dx,
    lx=lx,Lx=Lx,Tx=Tx,ex=ex
    )
  }

