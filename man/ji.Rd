\name{ji}
\alias{ji}
\title{ji Jaccard Index}
\description{The Jaccard similarity (overlap) index measures the size of the intersection of two populations relative to size of their union. It varies between zero (no overlap) and one (perfect overlap). The Jaccard Index is closely related the Sorensen (implemented in function \code{\link{li}}) and the Dice indices which are widely used in both the ecological and immunological literature (see, Rempala and Seweryn 2013).
}
\usage{ji(x, CI = 0.95, resample = 100, graph = FALSE, csv_output = FALSE, 
PlugIn = FALSE, size = 1, saveBootstrap = FALSE)}
\arguments{
\item{x}{a matrix containing input populations} 
\item{CI}{Confidence Interval default = 0.95, range (0, 1)} 
\item{resample}{set number of repetitions, default = 100} 
\item{graph}{default = FALSE, plot the results of hierarchical clustering of pairwise analysis of Jaccard Index; \code{graph = 'fileName'} user-defined output file name}   
\item{csv_output}{save the result of the analysis as .CSV file, default = FALSE; \code{csv_output = 'fileName'} user-defined output file name } 
\item{PlugIn}{standard plug-in estimator, default = FALSE} 
\item{size}{resampled fraction of the population, default = 1 (actual size of populations). The value should not be smaller than 10\% of population (size = 0.1)}
\item{saveBootstrap}{Saves bootstrap result to a file. Use saveBootstrap = TRUE to save bootstrap results to a Bootstrap folder in current directory; saveBootstrap = 'FolderName' - saves bootstrap results to user-named folder}
}
\author{Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}}
\references{Rempala G.A., Seweryn M. (2013) Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68}
\examples{data(TCR.Data)
result <- ji(x, resample = 50)}

