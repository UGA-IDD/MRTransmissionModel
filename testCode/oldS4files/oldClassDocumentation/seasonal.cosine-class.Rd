\name{seasonal.cosine-class}
\Rdversion{1.1}
\docType{class}
\alias{seasonal.cosine-class}
\alias{get.seasonal.mult,numeric,seasonal.cosine-method}

\title{Class \code{"seasonal.cosine"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("seasonal.cosine", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{amplitude}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{offset}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{period}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{seasonal.force}"}, directly.
}
\section{Methods}{
  \describe{
    \item{get.seasonal.mult}{\code{signature(t.in.yr = "numeric", season.obj = "seasonal.cosine")}: ... }
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
showClass("seasonal.cosine")
}
\keyword{classes}
