\name{dp.ht}
\alias{dp.ht}
\title{dp.ht Diversity Profile with the Horvitz-Thompson Adjustment}
\description{Calculates diversity profile with the Horvitz-Thompson adjustment (DP-HT), as defined in Rempala and Seweryn (2013) using the Renyi entropy (Renyi 1961) as a diversity measure. The function calculates the Renyi entropy values for a given range of the Renyi index (the index should be greater than 0). When the index is less then one, the rare counts are up-weighted and when it is greater than one, the rare counts are down-weighted. Since the Renyi entropy is a non-increasing function of the index, the profile plot should be always non-increasing. For more information, see Rempala and Seweryn (2013).
}
\usage{dp.ht(x, alpha = seq(0.1, 2, 0.1), CI = 0.95, resample = 100, 
single_graph = FALSE, pooled_graph = FALSE, csv_output = FALSE, 
PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE)}
\arguments{
\item{x}{a matrix containing input populations} 
\item{alpha}{a vector containing alpha values, default = seq(0.1, 2, 0.1)}  
\item{CVG}{a vector containing alpha values multiplied by coverage; default = FALSE}
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
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}
}
\references{Rempala G.A., Seweryn M. (2013) Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68

Renyi P. (1961) On measures of information and entropy. In: Proceedings of the 4th Berkeley symposium
on mathematics, statistics and probability, pp 547-61

Tothmeresz B. (1995) Comparison of different methods for diversity ordering. J Veget Sci 6:283-90

}
\examples{data(TCR.Data)
result <- dp.ht(x, PlugIn = TRUE)}

