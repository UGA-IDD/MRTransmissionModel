\name{experiment.updatedemog.vaccinationchange.vaccinationlimitations-class}
\Rdversion{1.1}
\docType{class}
\alias{experiment.updatedemog.vaccinationchange.vaccinationlimitations-class}
\alias{run,experiment.updatedemog.vaccinationchange.vaccinationlimitations-method}

\title{Class \code{"experiment.updatedemog.vaccinationchange.vaccinationlimitations"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("experiment.updatedemog.vaccinationchange.vaccinationlimitations", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{MR1MR2correlation}:}{Object of class \code{"logical"} ~~ }
    \item{\code{MR1SIAcorrelation}:}{Object of class \code{"logical"} ~~ }
    \item{\code{MR2SIAcorrelation}:}{Object of class \code{"logical"} ~~ }
    \item{\code{SIAinefficient}:}{Object of class \code{"logical"} ~~ }
    \item{\code{SIAinacc}:}{Object of class \code{"logical"} ~~ }
    \item{\code{prop.inacc}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{time.specific.MR1cov}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.MR2cov}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.SIAcov}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.min.age.MR1}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.max.age.MR1}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.min.age.MR2}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.max.age.MR2}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.min.age.SIA}:}{Object of class \code{"vector"} ~~ }
    \item{\code{time.specific.max.age.SIA}:}{Object of class \code{"vector"} ~~ }
    \item{\code{obj.vcdf.MR1}:}{Object of class \code{"vaccine.cdf.byage"} ~~ }
    \item{\code{obj.vcdf.MR2}:}{Object of class \code{"vaccine.cdf.byage"} ~~ }
    \item{\code{obj.prob.vsucc}:}{Object of class \code{"prob.vsucc.byage"} ~~ }
    \item{\code{sia.timing.in.year}:}{Object of class \code{"vector"} ~~ }
    \item{\code{trans}:}{Object of class \code{"ID.transition.SIR.vac"} ~~ }
    \item{\code{surv.each.timestep}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pop.rescale.each.timestep}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{maternal.obj}:}{Object of class \code{"maternal.exp.decay"} ~~ }
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
Class \code{"\linkS4class{experiment.updatedemog.vaccinationchange}"}, directly.
Class \code{"\linkS4class{experiment.updatedemog}"}, by class "experiment.updatedemog.vaccinationchange", distance 2.
Class \code{"\linkS4class{experiment}"}, by class "experiment.updatedemog.vaccinationchange", distance 3.
}
\section{Methods}{
  \describe{
    \item{run}{\code{signature(exper = "experiment.updatedemog.vaccinationchange.vaccinationlimitations")}: ... }
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
showClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations")
}
\keyword{classes}
