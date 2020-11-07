# Comparative-tour-through-simulation-algorithms-for-max-stable-processes

This is the repository for the code related to the paper

M. Oesting & K. Strokorb (2020+), A comparative tour through the simulation algorithms for max-stable processes, *arXiv preprint* [1809.09042](https://arxiv.org/abs/1809.09042).

The repository consists of two folders **Section2_DataExample** and **Section6_NumericalResults** containing the code related to the corresponding sections of the paper.

## Section2_DataExample

The folder **Section2_DataExample** contains the following six files:
- *DataExample.R*
- *DataExample_Advanced.R*
- *Algorithms.R*
- *AuxiliaryFunctions.R*
- *Section2.RData*
- *inlandBRsimu_0001.RData*

*DataExample.R* allows to retrace all estimation and simulation steps in Section 2 ``Data Example'' of the manuscript. The actual simulation procedure used therein is part of *Algorithms.R* and sourced from *DataExample.R*. Please find a detailed description of the functions in *Algorithms.R* in the folder **Section6_NumericalResults** and the *README_algorithms.txt* file therein.

*DataExample_Advanced.R* is an extension of *DataExample.R* that generates several intermediate illustrative plots along the way. It relies on *AuxiliaryFunctions.R* and loads more `R` packages related to the illustration of spatial data.

Both implementations (*DataExample.R* and its extension *DataExample_Advanced.R*) load the original temperature data under investigation and some auxiliary quantities from *Section2.RData*.

Since the simulation of 60000 samples from the fitted max-stable process on the considered inland grid takes several days on a standard laptop, we have exemplarily included our first 600 simulations in the file *inlandBRsimu_0001.RData*.
Thereby, the reader will be able to retrace how to obtain answers for questions Q1 and Q2 in the article as illustrated here based on these (insufficient) 600 samples.

## Section6_NumericalResults

The folder **Section6_NumericalResults** contains the following four files:
- *Algorithms.R*
- *README_algorithms.txt*
- *Code_SimulationStudy.R*
- *SimulationStudy_Param.RData*

The three generic simulation procedures used to obtain the numerical results are included in the file *Algorithms.R*. A detailed description of the functions therein is given in *README_algorithms.txt*.

*Code_SimulationStudy.R* contains sample code to recover the results of the numerical experiments for both Brown-Resnick and extremal-t processes. The input parameters for the simulation procedures in the different settings are given in *SimulationStudy_Param.RData*.

## Data

The file *Section2.RData* contains daily maximum temperatures from 1990 to 2019 that were measured at 18 inland stations in the Netherlands and are freely available from [knmi.nl](http://projects.knmi.nl/klimatologie/daggegevens/selectie.cgi).

## Third-Party Code

Some auxiliary functions for the simulation are taken from and parts of simulation algorithms themselves are based on the supplementary material of

C. Dombry, S. Engelke & M. Oesting (2016), Exact simulation of max-stable processes, *Biometrika* 103(2), pp.303-317 

available from [academic.oup.com](https://doi.org/10.1093/biomet/asw008).
