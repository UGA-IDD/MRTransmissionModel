\name{seasonal.periodic.bspline-class}
\Rdversion{1.1}
\docType{class}
\alias{seasonal.periodic.bspline-class}
\alias{get.seasonal.mult,numeric,seasonal.periodic.bspline-method}

\title{Class \code{"seasonal.periodic.bspline"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("seasonal.periodic.bspline", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{parameters}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nbasis}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{degree}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{period}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{seasonal.force}"}, directly.
}
\section{Methods}{
  \describe{
    \item{get.seasonal.mult}{\code{signature(t.in.yr = "numeric", season.obj = "seasonal.periodic.bspline")}: ... }
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
showClass("seasonal.periodic.bspline")
}
\keyword{classes}
