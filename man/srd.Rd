\name{srd}
\alias{srd}
\title{srd Symmetrized Renyi's Divergence}
\description{The symmetrized Renyi divergence (RD) is a measure of similarity between two discrete probability distributions which is non negative and symmetric. For details, see the description of function \code{\link{rd}} or Rempala and Seweryn (2013).}
\usage{srd(x, alpha = 0.5, CI = 0.95, resample = 100, graph = FALSE, csv_output = FALSE, 
PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE)}
\arguments{
\item{x}{a matrix containing input populations} 
\item{alpha}{Renyi's Divergence index of order \emph{alpha} must be between 0 and 1, default = 0.5}  
\item{CVG}{Renyi's Divergence index of order \emph{alpha} = coverage. If CVG = TRUE argument \emph{alpha} is ignored; default = FALSE}
\item{CI}{Confidence Interval default = 0.95, range (0, 1)} 
\item{resample}{number of repetitions, default = 100} 
\item{graph}{default = FALSE, plots the results of hierarchical clustering of pairwise analysis of Renyi's Divergence; \code{graph = 'fileName'} user-defined output file name}   
\item{csv_output}{save the result of the analysis as .CSV file, default = FALSE; \code{csv_output = 'fileName'} user-defined output file name } 
\item{PlugIn}{standard plug-in estimator, default = FALSE} 
\item{size}{resampled fraction of the population, default = 1 (actual size of populations). The value should not be smaller than 10\% of population (size = 0.1)}
\item{saveBootstrap}{Saves bootstrap result to a file. Use saveBootstrap = TRUE to save bootstrap results to a Bootstrap folder in current directory; saveBootstrap = 'FolderName' - saves bootstrap results to user-named folder}
}
\author{Christoph Sadee, Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}}
\references{Rempala G.A., Seweryn M. (2013) Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68}
\examples{data(TCR.Data)
result <- srd(x, resample = 20, alpha=0.5)}

