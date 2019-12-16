\name{cvg}
\alias{cvg}
\title{cvg Coverage}
\description{Calculates the sample coverage estimate using the Good-Turing formula. The sample coverage is an estimate of the probability of pulling a new species in the next draw, given a set of past observations For more details on CVG see Good I.J. (1953).}
\usage{cvg(x)}
\arguments{
\item{x}{a vector containing input population} 
}
\author{Christoph Sadee, Maciej Pietrzak, Michal Seweryn, Cankun Wang, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}}
\references{
Good I.J. The population frequencies of species and the estimation of population parameters. (1953) Biometrika 40:237-64\cr

Rempala G.A., Seweryn M. Methods for diversity and overlap analysis in T-cell receptor populations. (2013) J Math Biol 67:1339-68
}
\examples{data(TCR.Data)
result <- cvg(x[,1])}

