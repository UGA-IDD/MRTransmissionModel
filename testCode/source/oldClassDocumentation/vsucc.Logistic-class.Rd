\name{vsucc.Logistic-class}
\Rdversion{1.1}
\docType{class}
\alias{vsucc.Logistic-class}
\alias{pvacsuccess,numeric,vsucc.Logistic-method}

\title{Class \code{"vsucc.Logistic"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("vsucc.Logistic", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{intercept}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{mo.eff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{full.efficacy}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{vsucc}"}, directly.
}
\section{Methods}{
  \describe{
    \item{pvacsuccess}{\code{signature(age = "numeric", vsucc.obj = "vsucc.Logistic")}: ... }
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
showClass("vsucc.Logistic")
}
\keyword{classes}
