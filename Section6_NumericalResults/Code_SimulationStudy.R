source("Algorithms.R")

## Load parameters for the simulation study
load("SimulationStudy_Param.RData")
## this file contains two data frames named "param.BR" and "param.ET"
## which consist of the parameter sets for the numerical experimanents
## on Brown-Resnock and extreml-t processes, respectively
 
#######################################
####### Brown-Resnick processes #######
##### (cf. Figure 5 in the paper) #####
#######################################

## each line of "param.BR" contains the parameters "alpha" and "variance"
## of the corresponding (semi-)variogram
##     gamma(h) = variance * |h|^alpha
## as well as the targeted error (tolerated error probability)
## "target.err" and the corresponding
## thresholds for the sum- and supnormalized spectral representation,
## rel.thresh.sumnormal and rel.thresh.supnormal, and the size of the subgrid
## for the extremal functions approach, rel.subgrid.

## for instance, simulation parameters for the semi-variogram 
##     gamma(h) = 2+|h|^0.6
## and targeted error 0.05 are contained in row 17 of the data frame;
## the results of our experiments for this seeting can be recovered via

coord <- seq(-1,1,0.004)
no.simu <- 50000
row.index <- 17

vario <- function(x) param.BR[row.index,"variance"]*abs(x)^param.BR[row.index,"alpha"]

## sum-normalized representation (Dieker-Mikosch)
res.sumnormal <- simu_sumnormal(no.simu=no.simu, coord=coord, vario=vario,
                                type="brownresnick", 
                                rel.thresh=param.BR[row.index,"rel.thresh.sumnormal"],
                                calc.err=TRUE)
cat("mean time for DM in terms of (Gaussian) spectral functions:",
    mean(res.sumnormal$Gauss.counter), "\n")
cat("error probability (for control)", mean(res.sumnormal$err), "\n")

## sup-normalized representation
res.supnormal <- simu_supnormal(no.simu=no.simu, coord=coord, vario=vario,
                                type="brownresnick", 
                                rel.thresh=param.BR[row.index,"rel.thresh.supnormal"],
                                calc.err=TRUE)
cat("mean time for SN in terms of (Gaussian) spectral functions:",
    mean(res.supnormal$Gauss.counter), "\n")
cat("error probability (for control)", mean(res.supnormal$err), "\n")

## extremal functions
res.extrfcts <- simu_extrfcts(no.simu=no.simu, coord=coord, vario=vario,
                              type="brownresnick", 
                              rel.subgrid=param.BR[row.index,"rel.subgrid"],
                              calc.err=TRUE)
cat("mean time for EF in terms of (Gaussian) spectral functions:",
    mean(res.extrfcts$Gauss.counter), "\n")
cat("error probability (for control)", mean(res.extrfcts$err), "\n")

#######################################
######## Extremel-t processes #########
##### (cf. Figure 7 in the paper) #####
#######################################

## each line of "param.ET" contains the parameters "dof" and "scale"
## of the corresponding extremal-t process with correlation
##     rho(h) =  exp(-|h|/scale)
## and nu = dof degrees of freedom
## as well as the targeted error (tolerated error probability)
## "target.err" and the corresponding
## thresholds for the sum- and sup-normalized spectral representation,
## rel.thresh.sumnormal and rel.thresh.supnormal, and the size of the subgrid
## for the extremal functions approach, rel.subgrid.

## for instance, simulation parameters for the correlation
##     rho(h) = exp(-|h/0.5)
## and nu = 2 degrees of freedom
## and targeted error 0.1 are contained in row 12 of the data frame;
## the results of our experiments for this seeting can be recovered via

coord <- seq(-1,1,0.004)
no.simu <- 50000
row.index <- 12

corr <- function(x) exp(- abs(x/param.ET[row.index,"scale"]))
dof  <- param.ET[row.index,"dof"]


## sum-normalized representation
res.sumnormal <- simu_sumnormal(no.simu=no.simu, coord=coord, corr=corr, dof=dof,
                                type="extremalt", 
                                rel.thresh=param.ET[row.index,"rel.thresh.sumnormal"],
                                calc.err=TRUE)
cat("mean time for DM:", mean(res.sumnormal$Gauss.counter), "\n")
cat("error probability (for control)", mean(res.sumnormal$err), "\n")

## sup-normalized representation
res.supnormal <- simu_supnormal(no.simu=no.simu, coord=coord, corr=corr, dof=dof,
                                type="extremalt", 
                                rel.thresh=param.ET[row.index,"rel.thresh.supnormal"],
                                calc.err=TRUE)
cat("mean time for SN in terms of (Gaussian) spectral functions:",
    mean(res.supnormal$Gauss.counter), "\n")
cat("error probability (for control)", mean(res.supnormal$err), "\n")

## extremal functions
res.extrfcts <- simu_extrfcts(no.simu=no.simu, coord=coord, corr=corr, dof=dof,
                              type="extremalt", 
                              rel.subgrid=param.ET[row.index,"rel.subgrid"],
                              calc.err=TRUE)
cat("mean time for EF in terms of (Gaussian) spectral functions:",
    mean(res.extrfcts$Gauss.counter), "\n")
cat("error probability (for control)", mean(res.extrfcts$err), "\n")