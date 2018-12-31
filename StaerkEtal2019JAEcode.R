# ============================================================================
# FILE NAME:    StaerkEtal2019JAEcode.R
# DATE CREATED: 30 Dec. 2018
# AUTHORS:      Fernando Colchero and Johanna Staerk
# DESCRIPTION:  Simulation of survival and fecundity trajectories and 
#               calculation of generation time proxies.
# NOTES: 
#  - Code used in Staerk et al (2019) J Appl Ecol.
#  - This code simulates 15 mortality and 15 fecundity trajectories. From the 
#    combinations of these trajectories 225 lifetables are constructed. The
#    fecundity values for all 225 simulated life tables are adjusted such that 
#    the resulting population growth rate could take three values,namely 
#    λ = 0.9, λ = 1, and λ = 1.1, resulting in 675 populations that were either 
#    declining, stationary or increasing. Seven approximations of generation 
#    time are calulated for the 675 populations.

# INDEX:
# ------
#     A) FUNCTIONS ..................................................... L  34
#          A.1) MORTALITY .............................................. L  38
#          A.2) FECUNDITY .............................................. L  75
#          A.3) FUNCTIONS TO HANDLE DEMOGRAPHIC RATES .................. L  96
#     B) CREATE AGE-SPECIFIC SURVIVAL AND FECUNDITY PROFILES ........... L 194
#          B.1) MORTALITY PARAMETERS ................................... L 196
#          B.2) FECUNDITY PARAMETERS ................................... L 243
#     C) CREATE MATRIX OF GENERATION TIME APPROXIMATIONS ............... L 278
#          C.1) MATRIX OF COMBINATIONS OF MORTALITY AND FECUNDITY ...... L 280
#          C.2) LEVELS OF LAMBDA ....................................... L 285
#          C.3) OUTPUT MATRICES ........................................ L 291
#          C.4) RUN LOOP TO CALCULATE GEN. TIME APPROXIMATIONS ......... L 306
# ==============================================================================

# ============
# A) FUNCTIONS 
# ============

# ---------------
# A.1) MORTALITY:
# ---------------
# MORTALITY OR HAZARDS RATE AT AGE x: The function can produce declining morta-
# lity with age for pre-adults, and constant, exponential (i.e. senescent), or 
# logistic mortality for adults. 
CalcMort <- function(x, theta) {
  # x = age
  # theta = vector of mortality parameters that define the shape of the
  # mortality trjectory
  theta[1] * exp(-theta[2] * x) + theta[3] + theta[4] * exp(theta[5] * x) / 
    (1 + theta[6] * theta[4] / theta[5] * (exp(theta[5] * x) - 1))
}

# CUMMULATIVE SURVIVAL FUNCTION:
CalcSurv <- function(x, theta) {
  # x = age
  # theta = vector of mortality parameters
  if (theta[6] != 0) {
    sx <- exp(theta[1] / theta[2] * (exp(-theta[2] * x) - 1) - theta[3] * x) * 
      (1 + theta[6] * theta[4] / theta[5] * 
         (exp(theta[5] * x) - 1))^(-1 / theta[6])
  } else {
    sx <- exp(theta[1] / theta[2] * (exp(-theta[2] * x) - 1) - theta[3] * x + 
                theta[4] / theta[5] * (1 - exp(theta[5] * x)))
  }
  return(sx)
}

# LIFE EXPECTANCY AT BIRTH: using a left-point approximation.
CalcLifeExp <- function(x, theta, Dx) {
  # x = age 
  # theta = vector of mortality parameters
  # Dx = width of the age interval
  sum(CalcSurv(x, theta) * Dx)
}

# -------------------------
# A.2) FECUNDITY:
# -------------------------
# FECUNDITY FUNCTION:
CalcFecund <- function(x, gamma, alpha = 1) {
  # x = age
  # gamma = vector of fecundity parameters
  # alpha = age at sexual maturity
  xx <- x[-c(1:alpha)] - alpha
  c(rep(0, alpha), gamma[1] * exp(gamma[2] * (xx - gamma[3])^2))
}

# FUNCTION TO STRETCH FECUNDITY FOR MOST OF THE LIFESPAN: 
CalcFecPars <- function(i, j) {
  fecPars <- c(1, 0, 0)
  fecPars[3] <- diff(durmat[j, c(3, 2)]) * mxMaxAge[i]
  x0 <- ifelse(i < 10, durmat[j, 2], durmat[j, "alpha"])
  fecPars[2] <- log(fecPars[1] * 0.01 / fecPars[1]) / (x0 - fecPars[3])^2
  return(fecPars)
}

# -------------------------------------------
# A.3) FUNCTIONS TO HANDLE DEMOGRAPHIC RATES:
# -------------------------------------------
# FUNCTION TO FILL-IN LESLIE MATRICES:
FillMatr <- function(p, f, n) {
  # p = vector of age specific survival
  # f = vector of age specific fecundity
  # n = number of ages
  idcol <- 1:n
  idrow <- c(2:n, n)
  idpx <- (idcol - 1) * n + idrow
  Ax <- matrix(0, n, n)
  Ax[1, ] <- f
  Ax[idpx] <- p
  return(Ax)
}

# FUNCTION TO EXTRACT AND RE-SCALE VITAL RATES: extracts vital rates, re-scales
# them, calculate corresponding constant rates, provide lambda values and 
# assymptotic age distributions (right eigenvector of age-specific 
# Leslie matrix). 
# Arguments:
#    - id: index for the combination of mortality and fecundity (depends on the
#          number of rows of the index matrix 'ij', see below)
#    - forceLam: can take NA or a numerical value for the target population 
#      growth rate. 
ExtractRates <- function(id, forceLam = NA) {
  i <- ij[id, 1]
  j <- ij[id, 2]
  x <- 0:ceiling(durmat[i, 2])
  px <- CalcSurv(x + 1, thepars[i, ]) / CalcSurv(x, thepars[i, ])
  alpha <- floor(durmat[i, 1] / 4)
  fecPars <- CalcFecPars(i, j) 
  mx <- CalcFecund(x, fecPars, alpha = alpha)
  nx <- length(px)
  Ax <- FillMatr(px, mx, nx)
  eAx <- eigen(Ax)
  lamDx <- Re(eAx$value[1])
  if (!is.na(forceLam)) {
    jj <- 0
    while(jj < 5) {
      jj <- jj + 1
      if (jj > 1 & lamDx > forceLam) {
        px <- px * 0.95
      }
      mx <- mx * (forceLam / lamDx)
      Ax <- FillMatr(px, mx, nx)
      eAx <- eigen(Ax)
      lamx <- Re(eAx$values[1])
      lcount <- 0
      lamold <- lamx
      nlamsame <- 0
      while(round(lamx, 4) != forceLam & lcount < 200 & nlamsame < 3) {
        lcount <- lcount + 1
        mx <- mx * forceLam / lamx
        Ax[1, ] <- mx
        eAx <- eigen(Ax)
        lamx <- Re(eAx$value[1])
        if (abs(lamold - lamx) < 0.0001) nlamsame <- nlamsame + 1
        lamold <- lamx
      }
    }
    lamDx <- lamx
  }
  ageDist <- abs(Re(eAx$vector[, 1]))
  ageDist <- ageDist / sum(ageDist)
  idJ <- 1:alpha
  idA <- (alpha + 1):nx
  ageDistJ <- ageDist[idJ] / sum(ageDist[idJ])
  ageDistA <- ageDist[idA] / sum(ageDist[idA])
  pc <- c(px[idJ], 
          rep(px[idA] %*% ageDistA, nx - alpha))
  mc <- c(rep(0, alpha), 
          rep(mx[idA] %*% ageDistA, nx - alpha))
  nc <- length(pc)
  Ac <- FillMatr(pc, mc, nc)
  eAc <- eigen(Ac)
  lamDc <- Re(eAc$value[1])
  return(list(px = px, mx = mx, pc = pc, mc = mc, lamx = lamDx, lamc = lamDc,
              AD = ageDist, ADJ = ageDistJ, ADA = ageDistA, x = x))
}

# GENERATION LENGTH FUNCTION: based on the 'period definition', namely, the 
# mean age of the parents of offspring produced in the current time period 
# once the population has reached the stable age distribution (Leslie 1966).
# However, depending on the values of p, f, and l, other proxies can be calcu-
# lated.
CalcGenLen <- function(x, p, f, l) {
  # x = age vector
  # p = vector of age-specific survival
  # f = vector of age-specific fecundity
  # l = assymptotic population growth rate, lambda.
  lx <- c(1, cumprod(p))
  lx <- lx[-length(lx)]
  Tx <- sum(x * l^(-x) * lx * f) /sum(l^(-x) * lx * f)
  return(Tx)
}

# ======================================================
# B) CREATE AGE-SPECIFIC SURVIVAL AND FECUNDITY PROFILES 
# ======================================================
# B.1) MORTALITY PARAMETERS:
# --------------------------
# Width of age interval
Dx <- 0.01

# Ages
x <- seq(0, 1000, Dx)

# Matrix of mortality parameters for 5 general mortality profiles: 
nmini <- 5
theParsIni <- matrix(0, nmini, 6)
theParsIni[1, ] <- c(0.05, 0.12, 0.075, 0.015, 0.015, 0)
theParsIni[2, ] <- c(0.1, 0.25, 0.005, 0.02, 0.08, 0)
theParsIni[3, ] <- c(0, 0.15, 0.025, 0.025, 0.09, 0)
theParsIni[4, ] <- c(0, 0.15, 0.05, 0.025, 0.1, 1)
theParsIni[5, ] <- c(0, 0.15, 0.075, 0.005, 0.2, 3)

# Matrix of 15 sets of mortality parameters based on those in theParsIni. Each
# row of parameters in theParsIni is scaled for low, mid, and high mortalities:
thepars <- theParsIni[0, ]
for (i in 1:nmini) {
  thepars <- rbind(thepars, theParsIni[i, ] * 2,
                 theParsIni[i, ], theParsIni[i, ] * 0.4)
}

# Number of sets of mortality parameters (i.e. nm = 15):
nm <- nrow(thepars)

# Matrix of life expectancy at birth and age at which only 0.001% of 
# individuals are still alive for the 15 mortality trajectories (used later 
# to stretch fecundity for most of the lifespan):
durmat <- matrix(0, nm, 3, 
                 dimnames = list(NULL, c("e0", "x:S(x)=0.001", "alpha")))
for (i in 1:nm) {
  durmat[i, 1] <- CalcLifeExp(x, thepars[i, ], Dx)
  durmat[i, 2] <- x[which(CalcSurv(x, thepars[i, ]) < 0.001)[1]]
  durmat[i, 3] <- floor(durmat[i, 1] / 2)
}

# Plot 15 mortality trajectories 
par(mfrow = c(5, 3), mar = c(3.1, 3.1, 1, 1))
for (i in 1:nm) {
  plot(x[x < durmat[i, 2]], CalcMort(x[x < durmat[i, 2]], thepars[i, ]), 
       type = 'l', xlab = "", ylab = "", xlim = c(0, max(durmat[, 2])))
  mtext(paste("e0 =", signif(durmat[i, 1], 2)), 3, line = -2)
}

# B.2) FECUNDITY PARAMETERS:
# --------------------------
# Matrix of baseline fecundity paramenters
f0 <- rep(1.2, 5)
f1 <- c(-0.01, -0.05, -0.05, -0.0025, -0.0005)
f2 <- c(-0.5, 4, 8, 25, 40)
gamParsIni <- cbind(f0, f1, f2)

# Number of initial sets of fecundity parameters:
nfini <- nrow(gamParsIni)

# Matrix of 15 sets of fecundity parameters based on those in gamParsIni. Each
# row of parameters in gamParsIni is scaled for low, mid, and high fecundities:
gamPars <- gamParsIni[0, ]
for (i in 1:nfini) {
  gamPars <- rbind(gamPars, gamParsIni[i, ] * c(3, 1.5, 1), 
                 gamParsIni[i, ], 
                 gamParsIni[i, ] * c(0.5, 0.5, 1))
}

# Number of sets of fecundity parameters (i.e. nf = 15):
nf <- nrow(gamPars)

# Vector of ages at the mode (i.e. maximum) of fecundity. Negative values imply
# that there is no mode (i.e. declining fecundity with age).
mxMaxAge <- rep(c(-0.2, 0.25, 0.5,0.75, 1.2), each = 3)

# Plot 15 fecundity trajectories
par(mfrow = c(5, 3), mar = c(3.1, 3.1, 1, 1))
for (i in 1:nf) {
  plot(x, CalcFecund(x, gamPars[i, ]), type = 'l', 
       xlab = "", ylab = "", xlim = c(0, 40))
}

# ---------------------------------------------------
# C) CREATE MATRIX OF GENERATION TIME APPROXIMATIONS 
# ---------------------------------------------------
# C.1) MATRIX OF COMBINATIONS OF MORTALITY AND FECUNDITY:
# Index of different combinations of survival and fecundity trajectories
ij <- expand.grid(i = 1:15, j = 1:15)
nij <- nrow(ij)

# C.2) LEVELS OF LAMBDA: Set three fixed levels of population growth rates; i.e.
# a) lambda = 0.9 -> declining pop; b) lambda = 1 -> stationary population;
# c) lambda = 1.1 -> increasing population.
lamlevs <- c(0.9, 1, 1.1)
nlam <- length(lamlevs)

# C.3) OUTPUT MATRICES:
# Total number of combinations tested:
niter <- nij * nlam

# Output matrix of generation time approximations:
Tmat <- matrix(NA, niter, 8, 
               dimnames = list(NULL, c("Tb", "Ts", "Tc", "Tcp",
                                       "Tcf", "Tsc", "Tz", "Tq")))

# Output matrix of ancillary data (i.e. lambda, reproductive lifespan, 
# age at first reproduction, age at last reproduction)
covDat <- matrix(NA, niter, 4,
                 dimnames = list(NULL, c("lambda", "reprLifesp", 
                                         "alpha", "omega")))

# C.4) RUN LOOP TO CALCULATE GEN. TIME APPROXIMATIONS:
# Progress bar 
progrBar <- txtProgressBar(min = 1, max = niter, style = 3)

# Initialize the general index:
ii <- 0

# Run Loops:
for (rcom in 1:nij) {
  for (ll in 1:nlam) {
    ii <- ii + 1
    xvr <- ExtractRates(rcom, forceLam = lamlevs[ll])
    nx <- length(xvr$px)
    px <- xvr$px # survival at age x
    pc <- xvr$pc # constant survival
    mx <- xvr$mx # fecundity at age x 
    mc <- xvr$mc # constant fecundity
    x <- xvr$x # age
    alpha <- x[min(which(mx > 0))] # age at first reproduction
    omega <- max(x) # age at last reproduction
    qc <- 1 - pc[alpha + 1] # constant mortality probability with age
    
    # Calculate lambda for constant survival and age-specific fecundity:
    Acp <- FillMatr(pc, mx, length(pc))
    eAcp <- eigen(Acp)
    lAcp <- Re(eAcp$values[1])
    
    # Calculate lambda for constant fecundity and age-specific survival:
    Acm <- FillMatr(px, mc, length(pc))
    eAcm <- eigen(Acm)
    lAcm <- Re(eAcm$values[1])
    
    # Calculate different generation time approximations:
    Tb <- CalcGenLen(x, px, mx, xvr$lamx)
    Tc <- CalcGenLen(x, pc, mc, xvr$lamc)
    Tcp <- CalcGenLen(x, pc, mx, lAcp)
    Tcm <- CalcGenLen(x, px, mc, lAcm)
    Ts  <- CalcGenLen(x, px, mx, 1)
    Tsc <- CalcGenLen(x, pc, mc, 1)
    
    # Store approximations in Tmat:
    Tmat[ii, ] <- c(Tb, Ts, Tc, Tcp, Tcm, Tsc, alpha + 0.28 * (omega - alpha), 
                    alpha + 1 / qc)
    
    # Store ancillary data:
    covDat[ii, ] <- c(xvr$lamx, omega-alpha, alpha, omega)
    
    # Print progress:
    setTxtProgressBar(progrBar, ii)
  }
}

# -----------
# END OF CODE
# -----------