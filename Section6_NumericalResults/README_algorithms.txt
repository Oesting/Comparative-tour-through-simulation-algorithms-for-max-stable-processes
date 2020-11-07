############################################################
The file "Algorithms.R" contains the three main functions
for simulation of max-stable processes (here: Brown-Resnick 
and extremal-t processes) according to
- the sum-normalized spectral representation 
  (Dieker & Mikosch 2015, Dombry, Engelke & Oesting 2016)
- the extremal fucntions approach
  (Dombry, Engelke & Oesting 2016)
- the sup-normalized spectral representation
  (Oesting, Schlather & Zhou 2018)
############################################################
Function Calls:
  
  simu_sumnormal(no.simu=1, coord, vario, corr, dof, type,
                 rel.thresh=1, calc.err=TRUE)

  simu_extrfcts(no.simu=1, coord, vario, corr, dof, type, 
                rel.subgrid=1, calc.err=TRUE)
  
  simu_supnormal(no.simu=1, coord, vario, corr, dof, type, 
                 rel.thresh=1, calc.err=TRUE)
############################################################
Arguments:
    
no.simu        integer; number of simulations to be performed
  
coord          vector or matrix consisting of the coordinates
               corresponding to the coordinates of the random
               vector; here, each row corresponds to one 
               component

vario          function; returns the value of the variogram 
               used for the Brown-Resnick model 

corr           function; returns the value of the correlation
               used for the extremal-t model

dof            positive number; the degrees of freedom used in
               the extremal-t model

type           either "brownresnick" or "extremalt" to specify
               whether the Brown-Resnick or the extrmel-t model
               is to be simulated

rel.thresh     a number between 0 and 1; the value of the threshold
               used for simulation viathe sum- and sup-normalized 
               representation relative to the threshold that is 
               required to exact simulation; consequently, a value 
               of 1 leads to exact simulation, while values smaller
               than 1 may lead to approximation errors.
              
rel.subgrid    a number between 0 and 1; the percentage of points in
               the subgrid, on which exact simulation is performed
               via the extremal functions approach; a value of 1 
               ensures exact simulation for the whole domain, while
               values smaller than 1 may lead to some approximation
               errors.

calc.err       logical; if TRUE, an exact simulation is computed
               and compared to the approximation obtained using
               rel.thresh and rel.subgrid, respectively.
#############################################################
Output:
  
The output is a list. If calc.err=FALSE, this list contains

res            A matrix including the (non-exact) simulation results.
               Each row corresponds to one simulation.
         
spec.counter   A vector consisting of the numbers of spectral processes
               simulated for each (non-exact) max-stable simulation.

Gauss.counter  A vector consisting of the numbers of Gaussian processes
               simulated for each (non-exact) max-stable simulation.
               Note that Gauss.counter and spec.counter are identical for
               simu_sumnormal and simu_extrfcts.
               
If calc.err=TRUE, the list also contains               

res.full            A matrix including the (theoretical) exact simulation 
                    results without errors. Each row corresponds to one 
                    simulation.

spec.counter.full   A vector consisting of the numbers of spectral processes
                    simulated for each exact max-stable simulation.

Gauss.counter.full  A vector consisting of the numbers of Gaussian processes
                    simulated for each exact max-stable simulation.
                    Note that Gauss.counter.full and spec.counter.full are
                    identical for simu_sumnormal and simu_extrfcts.