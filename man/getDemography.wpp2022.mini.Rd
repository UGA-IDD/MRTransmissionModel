% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDemography.wpp2022.mini.R
\name{getDemography.wpp2022.mini}
\alias{getDemography.wpp2022.mini}
\title{Function nested inside getDemography.wpp2022() to only pull pop and age.structure}
\usage{
getDemography.wpp2022.mini(
  uncode = NA,
  age.classes = c(1:101),
  if.missing.use.region = F
)
}
\arguments{
\item{uncode}{UN country code}

\item{age.classes}{age classes in years}

\item{if.missing.use.region}{logical; c(MHL, XK, TUV)}
}
\value{
demography for each country
}
\description{
Function nested inside getDemography.wpp2022() to only pull pop and age.structure
}
