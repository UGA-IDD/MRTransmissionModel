\name{ID.transition.SIR.vac-class}
\Rdversion{1.1}
\docType{class}
\alias{ID.transition.SIR.vac-class}
\alias{next.ID.state,ID.state.matrix,ID.transition.SIR.vac-method}

\title{Class \code{"ID.transition.SIR.vac"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ID.transition.SIR.vac", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{vac.per}:}{Object of class \code{"vacc.per.time.step.by.age"} ~~ }
    \item{\code{n.epi.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{epi.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{epi.class.label}:}{Object of class \code{"character"} ~~ }
    \item{\code{s.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{i.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{r.inds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.age.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{age.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{aging.rate}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{survival.rate}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{birth.rate}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{waifw}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{age.surv.matrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{introduction.rate}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{exponent}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{frequency.dep}:}{Object of class \code{"logical"} ~~ }
    \item{\code{is.stochastic}:}{Object of class \code{"logical"} ~~ }
    \item{\code{get.births}:}{Object of class \code{"function"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ID.transition.SIR}"}, directly.
}
\section{Methods}{
  \describe{
    \item{next.ID.state}{\code{signature(state = "ID.state.matrix", tran = "ID.transition.SIR.vac")}: ... }
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
showClass("ID.transition.SIR.vac")
}
\keyword{classes}
