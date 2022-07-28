\name{ID.state.matrix-class}
\Rdversion{1.1}
\docType{class}
\alias{ID.state.matrix-class}
\alias{next.ID.state,ID.state.matrix,ID.transition.MSIRV-method}
\alias{next.ID.state,ID.state.matrix,ID.transition.SIR-method}
\alias{next.ID.state,ID.state.matrix,ID.transition.SIR.vac-method}

\title{Class \code{"ID.state.matrix"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ID.state.matrix", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{n.epi.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{epi.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{epi.class.label}:}{Object of class \code{"character"} ~~ }
    \item{\code{n.age.class}:}{Object of class \code{"numeric"} ~~ }
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
    \item{next.ID.state}{\code{signature(state = "ID.state.matrix", tran = "ID.transition.MSIRV")}: ... }
    \item{next.ID.state}{\code{signature(state = "ID.state.matrix", tran = "ID.transition.SIR")}: ... }
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
showClass("ID.state.matrix")
}
\keyword{classes}
