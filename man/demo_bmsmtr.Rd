% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demo_utils.R
\name{demo_bmsmtr}
\alias{demo_bmsmtr}
\title{demo multi-source multi-task model training and interpretative plots}
\usage{
demo_bmsmtr(dY = 3)
}
\arguments{
\item{dY}{is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.}
}
\description{
\code{demo_bmsmtr} trains bmsmtr on random multi-source multi-task regression dataset.
}
\examples{

#for multi-source multi-task regression
demo_bmsmtr(dY = 3)

}
