\name{seasonal.age.specific-class}
\Rdversion{1.1}
\docType{class}
\alias{seasonal.age.specific-class}
\alias{get.seasonal.mult,numeric,seasonal.age.specific-method}

\title{Class \code{"seasonal.age.specific"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("seasonal.age.specific", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{amplitude}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{offset}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{period}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{age.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{chosen.age}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{seasonal.force}"}, directly.
}
\section{Methods}{
  \describe{
    \item{get.seasonal.mult}{\code{signature(t.in.yr = "numeric", season.obj = "seasonal.age.specific")}: ... }
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
showClass("seasonal.age.specific")
}
\keyword{classes}
