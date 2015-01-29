\name{x}
\alias{x}
\title{Example dataset for divo package}
\description{T-cell receptor repertoires sequenced using Ion Torrent technology. Dataset contains receptors found in four different organs, each with two functional populations (naive and regulatory (Treg)). Cells are isolated from colon (Col), peripheral lymph nodes (PLN), mesenteric lymph nodes (MLN) and thymus (Thym). TCR populations data are stored in a matrix (object named \code{x}). Each column of x contains sequenced counts of specific TCR variants in given organ population \code{x}
}

\references{Cebula A., Seweryn M., Rempala G.A., Pabla S.S., McIndoe R.A., Denning T.L., Bry L., Kraj P., Kisielow P., Ignatowicz L. (2013) Thymus-derived regulatory T cells contribute to tolerance to commensal microbiota. Nature 497:258-62}
\examples{data(TCR.Data)
head(x)}

