% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratefyR.R
\name{stratefyR}
\alias{stratefyR}
\title{CSR ecology strategies calculator}
\usage{
stratefyR(sp = NULL, la = NULL, lfw = NULL, ldw = NULL, data = NULL)
}
\arguments{
\item{sp}{species names or morphospecies codes}

\item{la}{leaf area (mm^2)}

\item{lfw}{leaf fresh weight (mg)}

\item{ldw}{leaf dry weight (mg)}

\item{data}{a data frame with species name, LA, LFW, LDW  values for each species.}
}
\value{
returns a data frame with C, S and R index strategies for each species.
}
\description{
This function...
}
\examples{

sp <- c("Molopospermum p.", "Kalmia p.", "Arabidopsis t.", "Pteridium a.")    #[A]
la <-  c(369615.7, 11.8,55.7,36061.2)                                          #[D]
lfw <- c(84009.6,4.5,12.3,6988.0)                                              #[E]
ldw <- c(21187.8,1.8,1.6,2480.0)                                               #[F]
traitsp <- data.frame(sp, la, lfw, ldw)
stratefyR(data = traitsp)    


}
\author{
Alexandre Adalardo de Oliveira \email{ecovirtualpackage@gmail.com}
}
\keyword{,}
\keyword{competition,}
\keyword{ecology}
\keyword{plants}
\keyword{ruderal,}
\keyword{strategies,}
\keyword{stress-tolerant}
