% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.prem.WAIFW.R
\name{get.prem.WAIFW}
\alias{get.prem.WAIFW}
\title{Make a WAIFW matrix based on Prem et al. 2021 using pakistan for afghanistan}
\usage{
get.prem.WAIFW(
  age.class.boundries = (1:90),
  uncode,
  other.contact.matrix = FALSE,
  bandwidth = c(3, 3),
  adjustment_start_year = FALSE,
  year = 1980
)
}
\arguments{
\item{age.class.boundries}{the upper age limit for each age class in YEARS}

\item{uncode}{country code}

\item{other.contact.matrix}{xxx}

\item{bandwidth}{desired smooth bandwidth - default=c(3,3)}

\item{adjustment_start_year}{for 202110gavi_v3 given that Shaun's R0 estimates are adjusted to 1980 pop structure}

\item{year}{xxx}
}
\value{
a WAIFW matrix based on the Polymod results from chosen location with row and colnames indicating age classes
}
\description{
Make a WAIFW matrix based on Prem et al. 2021 using pakistan for afghanistan
}
