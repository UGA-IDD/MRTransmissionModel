% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.births.fun.R
\name{get.births.fun}
\alias{get.births.fun}
\title{This function is an argument of the experiment builder to ensure fertility from susceptible mothers is susceptible}
\usage{
get.births.fun(
  state,
  tran,
  fert.curve = c(0, 36.3, 210.6, 164.8, 68.3, 27.5, 10.1, 2.8, 0),
  lower.age.boundary = c(0, 15, 20, 25, 30, 35, 40, 45, 50, 101)
)
}
\arguments{
\item{state}{state vector at this stage in simulation, length n.epi.class*n.age.class}

\item{tran}{transition matrix at this stage in simulation}

\item{fert.curve}{age specific fertility rate used to establish fraction susceptible mothers}

\item{lower.age.boundary}{corresponding lower age boundaries for the fert.curve}
}
\value{
vector of births of length state
}
\description{
This function is an argument of the experiment builder to ensure fertility from susceptible mothers is susceptible
}
