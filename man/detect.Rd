\name{detect}
\alias{detect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detect cancer classfication result}
\description{This function reads a patient's SNV patterns and return the classfication result}
\usage{
detect(x , cancer , k = 6 , ignoreSame = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{MAF file}
  \item{cancer}{The list of cancers that processed in getMatrix.R}
  \item{k}{The k of KNN}
  \item{ignoreSame}{Dicedes whether ignore the exactly same matrices in all six dimensions}
}
\details{
nothing
}
\value{
\item{which(classfication == max(classfication),arr.ind=TRUE)    }{classfication result}

{same    }{TRUE if there is a same sample}

}
\references{
nothing}
\author{
Bo Li}
\note{
nothing
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
nothing
}
\examples{
ki = system.file("extdata", "kirp.maf", package = "SNVCC")
pr = system.file("extdata", "prad.maf", package = "SNVCC")
kirp = getMatrix(ki)
prad = getMatrix(pr)
cancer = list(kirp,prad)
x = list(cancer[[1]][[1]][1,],cancer[[1]][[2]][1,],cancer[[1]][[3]][1,],cancer[[1]][[4]][1,],cancer[[1]][[5]][1,],cancer[[1]][[6]][1,])
detect(x , cancer , 6 , TRUE)
}

\keyword{ SNV }% use one of  RShowDoc("KEYWORDS")
\keyword{ Classfication }% __ONLY ONE__ keyword per line
