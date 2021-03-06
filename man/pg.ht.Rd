\name{pg.ht}
\alias{pg.ht}
\title{pg.ht Power-Geometric Index with the Horvitz-Thompson Correction}
\description{The Horvitz-Thompson corrected version of the Power Geometric (PG) index (see help for \code{\link{pg}}). The PG index is a generalization of the Morisita-Horn index as well as the Bhattacharyya's coefficient. It quantifies overlap as cosine of an angle between two exponentially normalized population vectors. For further details and definitions, see Rempala and Seweryn (2013).
}

\usage{pg.ht(x, alpha = 1, beta=alpha, CI = 0.95, resample = 100, graph = FALSE, 
csv_output = FALSE, PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE)}
\arguments{
\item{x}{a matrix containing input populations} 
\item{alpha}{PG of order \emph{alpha} < 1 puts more weight on the rare species and the I Index of order \emph{alpha} > 1 puts more weight on the abundant ones for first population, default = 1}  
\item{beta}{PG of order \emph{beta} < 1 puts more weight on the rare species and the I Index of order \emph{beta} > 1 puts more weight on the abundant ones for second population, default = alpha} 
\item{CVG}{PG of order \emph{alpha} or \emph{beta} = coverage. If CVG = TRUE argument \emph{alpha} is ignored; default = FALSE}
\item{CI}{Confidence Interval default = 0.95, range (0, 1)} 
\item{resample}{set number of repetitions, default = 100} 
\item{graph}{default = FALSE, plot the results of hierarchical clustering of pairwise analysis of Power-Geometric Index, \code{graph = 'fileName'} user-defined output file name}   
\item{csv_output}{save the result of the analysis as .CSV file, default = FALSE; \code{csv_output = 'fileName'} user-defined output file name } 
\item{PlugIn}{standard plug-in estimator, default = FALSE} 
\item{size}{resampled fraction of the population, default = 1 (actual size of populations). The value should not be smaller than 10\% of population (size = 0.1)}
\item{saveBootstrap}{Saves bootstrap result to a file. Use saveBootstrap = TRUE to save bootstrap results to a Bootstrap folder in current directory; saveBootstrap = 'FolderName' - saves bootstrap results to user-named folder}
}
\author{Christoph Sadee, Maciej Pietrzak, Michal Seweryn, Cankun Wang, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}}
\references{Rempala G.A., Seweryn M. (2013) Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68}
\examples{data(TCR.Data)
result <- pg.ht(x, PlugIn = TRUE)}

