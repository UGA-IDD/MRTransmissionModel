% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/space.wrapper.surv.prob.over.age.time.R
\name{space.wrapper.surv.prob.over.age.time}
\alias{space.wrapper.surv.prob.over.age.time}
\title{Wrapper function to add space to create.surv.prob.over.age.time()}
\usage{
space.wrapper.surv.prob.over.age.time(
  age.classes,
  generation.time,
  space.nMx = NULL,
  years.interpolate = seq(1980, (1980 + 121), 1),
  check = FALSE
)
}
\arguments{
\item{age.classes}{vector; a set of age classes for which tran is being built}

\item{generation.time}{numeric; generation time}

\item{space.nMx}{space.nMx object; age and space specific death rates over time}

\item{years.interpolate}{vector; years to interpolate age specific death rates}

\item{check}{boolean; whether or not to plot the death rates over age}
}
\value{
returns the age and spatial profile of survivorship (in units of the generation time) (rows) by time (cols)
}
\description{
Wrapper function to add space to create.surv.prob.over.age.time()
}
