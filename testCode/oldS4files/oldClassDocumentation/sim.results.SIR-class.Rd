\name{sim.results.SIR-class}
\Rdversion{1.1}
\docType{class}
\alias{sim.results.SIR-class}
\alias{getAverageInfectionAge,sim.results.SIR-method}
\alias{getCompartmentTotals,sim.results.SIR-method}
\alias{plot,sim.results.SIR-method}
\alias{plot,sim.results.SIR,ANY-method}

\title{Class \code{"sim.results.SIR"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sim.results.SIR", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{s.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{i.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{r.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{age.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{mMatrix}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{OptionalMatrix}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{OptionalNumericOrMatrix}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{MatrixorFunction}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{OptionalNumericOrMatrixOrCall}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{getAverageInfectionAge}{\code{signature(sim.res = "sim.results.SIR")}: ... }
    \item{getCompartmentTotals}{\code{signature(sim.res = "sim.results.SIR")}: ... }
    \item{plot}{\code{signature(x = "sim.results.SIR")}: ... }
    \item{plot}{\code{signature(x = "sim.results.SIR", NA = "ANY")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("sim.results.SIR")
}
\keyword{classes}
