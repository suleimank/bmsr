% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demo_utils.R
\name{demo_internal}
\alias{demo_internal}
\title{demo multi-source/multi-task model training and interpretative plots}
\usage{
demo_internal(file = "bmsr.stan", dY = 1)
}
\arguments{
\item{file}{is the stan file containing the stan code.}

\item{dY}{is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.}
}
\description{
\code{demo_bmsmtr} trains bmsr and bmsmtr on random multi-source/multi-task regression dataset.
}
\examples{

#for multi-source multi-task regression
demo_internal(file = "bmsr.stan",dY = 1)

#for multi-source multi-task regression
demo_internal(file = "bmsmtr.stan",dY = 3)

}
\keyword{internal}
