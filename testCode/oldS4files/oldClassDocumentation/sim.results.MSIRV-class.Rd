\name{sim.results.MSIRV-class}
\Rdversion{1.1}
\docType{class}
\alias{sim.results.MSIRV-class}
\alias{getCompartmentTotals,sim.results.MSIRV-method}
\alias{getCumRiskTrimester,sim.results.MSIRV,numeric-method}
\alias{plot,sim.results.MSIRV-method}
\alias{plot,sim.results.MSIRV,ANY-method}
\alias{toProportions,sim.results.MSIRV-method}

\title{Class \code{"sim.results.MSIRV"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sim.results.MSIRV", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{m.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{v.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{routine.intro}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sia.times}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{s.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{i.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{r.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{age.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sim.results.SIR}"}, directly.
Class \code{"\linkS4class{matrix}"}, by class "sim.results.SIR", distance 2.
Class \code{"\linkS4class{array}"}, by class "sim.results.SIR", distance 3.
Class \code{"\linkS4class{mMatrix}"}, by class "sim.results.SIR", distance 3.
Class \code{"\linkS4class{OptionalMatrix}"}, by class "sim.results.SIR", distance 3.
Class \code{"\linkS4class{OptionalNumericOrMatrix}"}, by class "sim.results.SIR", distance 3.
Class \code{"\linkS4class{MatrixorFunction}"}, by class "sim.results.SIR", distance 3.
Class \code{"\linkS4class{structure}"}, by class "sim.results.SIR", distance 4.
Class \code{"\linkS4class{OptionalNumericOrMatrixOrCall}"}, by class "sim.results.SIR", distance 4.
Class \code{"\linkS4class{vector}"}, by class "sim.results.SIR", distance 5, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{getCompartmentTotals}{\code{signature(sim.res = "sim.results.MSIRV")}: ... }
    \item{getCumRiskTrimester}{\code{signature(sim.res = "sim.results.MSIRV", nbiweeks = "numeric")}: ... }
    \item{plot}{\code{signature(x = "sim.results.MSIRV")}: ... }
    \item{plot}{\code{signature(x = "sim.results.MSIRV", NA = "ANY")}: ... }
    \item{toProportions}{\code{signature(ob = "sim.results.MSIRV")}: ... }
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
showClass("sim.results.MSIRV")
}
\keyword{classes}
