% Hand-maintained (not roxygen2-generated: there is no corresponding
% documented object in R/ for an Rcpp Modules reference class, so roxygen2
% has nothing to tag this with). Edit this file directly.
\name{C++Object-class}
\docType{class}
\alias{C++Object-class}
\alias{C++Object}
\title{Rcpp's generic C++ object base class}
\description{
\code{C++Object} is the S4 reference class that Rcpp Modules generates as the
common parent of every Rcpp Modules reference class, including this
package's own \code{\link[=Rcpp_neuron-class]{Rcpp_neuron}}. It is generated
by Rcpp itself, not by neuronsDG, and carries no behavior or documentation
specific to this package; it is documented here only because every exported
S4 class, including inherited ones that happen to be reachable this way,
needs a documentation entry for CRAN's checks. There is no reason for a user
of neuronsDG to reference \code{C++Object} directly.
}
\keyword{internal}
