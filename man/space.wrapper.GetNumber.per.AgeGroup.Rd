% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/space.wrapper.GetNumber.per.AgeGroup.R
\name{space.wrapper.GetNumber.per.AgeGroup}
\alias{space.wrapper.GetNumber.per.AgeGroup}
\title{Spatial wrapper Function to get number of individuals per age group and subpopulation}
\usage{
space.wrapper.GetNumber.per.AgeGroup(trans, state)
}
\arguments{
\item{trans}{transition object}

\item{state}{vector of state of the population at one time point}
}
\value{
matrix, size n.subpops (rows) by n.ages.class (columns) with the number of individuals in each age and space
}
\description{
Spatial wrapper Function to get number of individuals per age group and subpopulation
}
