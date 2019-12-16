\name{ens.ht}
\alias{ens.ht}
\title{ens.ht Effective Number of Species with the Horvitz-Thompson Correction}
\description{Calculates diversity profile (DP) using the effective number of species (ENS) based on inverting the Renyi entropy with the Horvitz-Thompson correction. For any monotone diversity index (see, e.g., Rempala and Seweryn 2013) the ENS is defined as the size of a uniform population with the same index value as the current population. The ENS may be considered as a measure of population diversity expressed in the units of species counts. The ENS profile is calculated against the Renyi entropy index, which allows for a direct comparison with the diversity profile (as in \code{\link{dp}}). The ENS without the Horvitz-Thompson correction is available as function \code{\link{ens}}. For more details on ENS see Rempala and Seweryn (2013) or Jost (2006).
}
\usage{ens.ht(x, alpha = seq(0.1, 2, 0.1), CI = 0.95, resample = 100, 
single_graph = FALSE, pooled_graph = FALSE, csv_output = FALSE, 
PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE)}
\arguments{
\item{x}{a matrix containing input populations} 
\item{alpha}{a vector containing alpha values, default = seq(0.1, 2, 0.1)}  
\item{CVG}{a list containing alpha values multiplied by coverage; default = FALSE}
\item{CI}{Confidence Interval default = 0.95, range (0, 1)} 
\item{resample}{set number of repetitions, default = 100} 
\item{single_graph}{default = FALSE, plot of the Diversity Profile for each population; \code{single_graph = 'fileName'} user-defined output file name}   
\item{pooled_graph}{default = FALSE, plot of the Diversity Profile for all populations; \code{pooled_graph = 'fileName'} user-defined output file name} 
\item{csv_output}{save the result of the analysis as .CSV file, default = FALSE; \code{csv_output = 'fileName'} user-defined output file name } 
\item{PlugIn}{standard plug-in estimator, default = FALSE} 
\item{size}{resampled fraction of the population, default = 1 (actual size of populations). The value should not be smaller than 10\% of population (size = 0.1)}
\item{saveBootstrap}{Saves bootstrap result to a file. Use saveBootstrap = TRUE to save bootstrap results to a Bootstrap folder in current directory; saveBootstrap = 'FolderName' - saves bootstrap results to user-named folder}
}
\author{Christoph Sadee, Maciej Pietrzak, Michal Seweryn, Cankun Wang, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}}
\references{
Jost L. (2006) Entropy and diversity. Oikos 113:363-75

Rempala G.A., Seweryn M. Methods for diversity and overlap analysis in T-cell receptor populations. (2013) J Math Biol 67:1339-68}
\examples{data(TCR.Data)
result <- ens.ht(x, PlugIn = TRUE)}

