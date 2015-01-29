\name{i.in}
\alias{i.in}
\title{i.in Information Index (I index) for 2-Way Table}
\description{The I-index is a measure of overlap in two way tables based on the generalized mutual information statistic. The I-index measures dependence among columns of two-way tables, taking values between 0 and 1. It returns a value of zero when the table columns form an orthogonal system and a value of one when the table columns rank is one. The value of the parameter alpha is related to the structure of dependence, as described in Rempala and Seweryn (2013).
}
\usage{i.in(x, alpha = 1, CI = 0.95, resample = 100, PlugIn = FALSE, size = 1, CVG = FALSE, 
saveBootstrap = FALSE)}
\arguments{
\item{x}{a matrix containing input populations} 
\item{alpha}{I index of order \emph{alpha}, must be between 0 and 1, default = 0.5}  
\item{CVG}{I index of order \emph{alpha} = coverage. If CVG = TRUE argument \emph{alpha} is ignored; default = FALSE}
\item{CI}{Confidence Interval default = 0.95, range (0, 1)} 
\item{resample}{set number of repetitions, default = 100} 
\item{PlugIn}{standard plug-in estimator, default = FALSE} 
\item{size}{resampled fraction of the population, default = 1 (actual size of populations). The value should not be smaller than 10\% of population (size = 0.1)}
\item{saveBootstrap}{Saves bootstrap result to a file. Use saveBootstrap = TRUE to save bootstrap results to a Bootstrap folder in current directory; saveBootstrap = 'FolderName' - saves bootstrap results to user-named folder}
}
\author{Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}}
\references{Rempala G.A., Seweryn M. (2013)  Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68}
\examples{data(TCR.Data)
result <- i.in(x, resample = 50)}

