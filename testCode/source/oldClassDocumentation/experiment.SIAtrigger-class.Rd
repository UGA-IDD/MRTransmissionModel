\name{experiment.SIAtrigger-class}
\Rdversion{1.1}
\docType{class}
\alias{experiment.SIAtrigger-class}

\title{Class \code{"experiment.SIAtrigger"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("experiment.SIAtrigger", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{trans}:}{Object of class \code{"ID.transition.SIR.vac"} ~~ }
    \item{\code{sia.obj}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{trigger.list}:}{Object of class \code{"list"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{R0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{state.t0}:}{Object of class \code{"ID.state.matrix"} ~~ }
    \item{\code{t.min}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t.max}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t0.doy}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{step.size}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{season.obj}:}{Object of class \code{"seasonal.force"} ~~ }
    \item{\code{births.per.1000.each.timestep}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{description}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{experiment}"}, directly.
}
\section{Methods}{
No methods defined with class "experiment.SIAtrigger" in the signature.
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
showClass("experiment.SIAtrigger")
}
\keyword{classes}
