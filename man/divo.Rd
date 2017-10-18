\name{divo}
\alias{divo}
\title{Tools for Analysis of Diversity and Similarity in Biological Systems}
\description{The package divo implements various algorithms for empirical analysis of diversity and similarity (overlap) in biological or ecological systems. The implemented indices of diversity and overlap are based both on the information-theoretic and geometric considerations. The indices have the capacity to naturally up-weight or down-weight rare and abundant population species counts, by applying the Good-Turing sample coverage correction. The functional version of a diversity index, the so-called diversity profile, is also implemented along with the diversity and overlap indices inversions known as the effective numbers of species (ENS).\cr
\strong{For examples and detailed information on specific functions, see their manual pages:}

\tabular{rlll}{
\tab\code{\link{cvg}}\tab    Coverage\cr
\tab\code{\link{dp}}\tab    Diversity Profile\cr
\tab\code{\link{dp.ht}}\tab    Diversity Profile with the Horvitz-Thompson Correction\cr
\tab\code{\link{ens}}\tab   Effective Number of Species\cr
\tab\code{\link{ens.ht}}\tab   Effective Number of Species with the Horvitz-Thompson Correction\cr
\tab\code{\link{i.in}}\tab    Information Index (I-index) for 2-Way Table\cr
\tab\code{\link{i.inp}}\tab    Information Index (I-index) for 2-Way, 2 Column Table\cr
\tab\code{\link{ji}}\tab    Jaccard Index\cr
\tab\code{\link{li}}\tab    Sorensen Index\cr
\tab\code{\link{mh}}\tab    Morisita-Horn Index\cr
\tab\code{\link{pg}}\tab    Power-Geometric Index\cr
\tab\code{\link{pg.ht}}\tab Power-Geometric Index with the Horvitz-Thompson Correction\cr
\tab\code{\link{rd}}\tab    Renyi's Divergence\cr
\tab\code{\link{srd}}\tab   Symmetrized Renyi's Divergence\cr
}
}

\references{
Rempala G.A., Seweryn M. (2013) Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68
}
\author{
Christoph Sadee, Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala\cr
Maintainer: Maciej Pietrzak \email{pietrzak.20@osu.edu}
}

\keyword{Contingency tables}
\keyword{Antigen receptors}
\keyword{Richness and diversity estimation}
\keyword{Renyi's entropy}
\keyword{Renyi's divergence}

