\name{getMatrix}
\alias{getMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read MAF files
}
\description{
This function reads MAF files and extracts six mutli-dimensional SNV matrices
}
\usage{
getMatrix(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
MAF file
}
}
\details{
nothing
}
\value{

{nu1 }{one-dimensional feature matrix}

{nu2_l }{two-dimensional left-adjacent feature matrix}

{nu2_r }{two-dimensional right-adjacent feature matrix}

{nu3_l }{three-dimensional left-adjacent feature matrix}

{nu3_m }{three-dimensional both sides adjacent feature matrix}

{nu3_r }{three-dimensional right-adjacent feature matrix}
{id }{id of MAF file}

}
\references{
nothing
}
\author{
Bo Li
}
\note{
nothing
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
nothing
}
\examples{
ki = system.file("extdata", "kirp.maf", package = "SNVCC")
getMatrix(ki)
}

\keyword{ MAF }% use one of  RShowDoc("KEYWORDS")
\keyword{ Feature Matrix }% __ONLY ONE__ keyword per line
