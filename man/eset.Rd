\name{eset}
\alias{eset}
\alias{esetSpot}
\alias{cloneanno}
\docType{data}
\title{renal cell cancer expression data}
\description{renal cell cancer expression data}
\usage{
data("eset")
data("esetSpot")
data("cloneanno")
}
\format{
\code{eset} is an object of type
\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.  
Samples (columns of the expression matrix) correspond to renal cell 
cancer samples from 74 patients with kidney tumors. Probes (rows of the
expression matrix) correspond to different cDNA clones.

\code{cloneanno} is a data frame with 4224 rows, one for each clone. Its
columns are:
\code{plate}, \code{SrcRow}, \code{SrcCol}, \code{imageid}, \code{AccNumber},
\code{spot1}, \code{spot2}, \code{description}, \code{vendor}

The arrays that were used to measure these data contained two features
for each reporter (clone), and \code{esetSpot} is the data \emph{before}
averaging over the two features per reporter. Spots 4225:8448 are the
duplicates of spots 1:4224.}

\source{Holger Sueltmann, DKFZ. Also available from 
ArrayExpress, Accession number E-DKFZ-1.

The data import fro Genepix files and Excel tables, and the
preprocessing is described in the files preproc.Rnw and preproc.pdf
in the \code{scripts} subdirectory of the package. Please refer to these for
details!
}

\references{
Gene expression in kidney cancer is associated with novel tumor subtypes, 
cytogenetic abnormalities and metastasis formation.
Holger Sueltmann, Anja von Heydebreck, Wolfgang Huber, Ruprecht Kuner, 
Andreas Buness, Markus Vogt, Bastian Gunawan, Martin Vingron, 
Laszlo Fuezesi, and Annemarie Poustka. Submitted.
}

\seealso{\code{\link{qua}}}
\keyword{datasets}
\examples{

data(eset)
eset
data(cloneanno)
cloneanno[1,]

## see also the vignette
}
