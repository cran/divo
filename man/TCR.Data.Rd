\name{TCR.Data}
\alias{TCR.Data}
\title{TCR.Data: Repertoires of Naive and Regulatory T-cell Populations}
\description{
T-cell receptor repertoires sequenced using Ion Torrent technology. 
Dataset contains receptors found in four different organs, each with two functional populations (naive and regulatory (Treg)). 
Cells are isolated from the colon (Col), peripheral lymph nodes (PLN), mesenteric lymph nodes (MLN), and thymus (Thym). 

TCR populations data are stored in a matrix (object named \code{x}). 
Each column of \code{x} contains sequenced counts of specific TCR variants in a given organ population.
}
\references{
Cebula A., Seweryn M., Rempala G.A., Pabla S.S., McIndoe R.A., Denning T.L., Bry L., Kraj P., Kisielow P., Ignatowicz L. (2013) 
Thymus-derived regulatory T cells contribute to tolerance to commensal microbiota. Nature 497:258-62.
}
\examples{
data(TCR.Data)
head(x)
}
