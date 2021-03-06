% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmsr.R
\name{predictSTAN}
\alias{predictSTAN}
\title{predict from a regression STAN model}
\usage{
predictSTAN(predFunction, out, xTest, nTest, yN = NULL)
}
\arguments{
\item{predFunction}{is a function for predicting y's using the outcome of stan run.}

\item{xTest}{is a matrix of test data for predicting the outcome. If NULL no prediction is made (default).}

\item{nTest}{is a vector of length S, containing the number of values in each source. Can contain zero's.}

\item{yN}{a list containing values used for normalizing the data: (default = NULL).
\itemize{
\item{cm}{vector of means with which the data is centered (0's if data is not centered),}
\item{cs}{vector of standard deviations with which the data is scaled (1's if data is not scaled).}
}}
}
\value{
yPred prediction vector of the stan model
}
\description{
\code{predictSTAN} predicts the output a stan regression model as defined by the parameter \code{predFunction}.
}
