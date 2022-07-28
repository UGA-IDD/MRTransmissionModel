\name{experiment-class}
\Rdversion{1.1}
\docType{class}
\alias{experiment-class}

\title{Class \code{"experiment"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("experiment", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{R0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{state.t0}:}{Object of class \code{"ID.state.matrix"} ~~ }
    \item{\code{t.min}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t.max}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t0.doy}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{step.size}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trans}:}{Object of class \code{"ID.transition.SIR"} ~~ }
    \item{\code{season.obj}:}{Object of class \code{"seasonal.force"} ~~ }
    \item{\code{births.per.1000.each.timestep}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{description}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "experiment" in the signature.
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
showClass("experiment")
}
\keyword{classes}
