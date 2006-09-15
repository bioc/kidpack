\name{qua}
\alias{qua}
\alias{spotanno}
\alias{hybanno}
\docType{data}
\title{renal cell cancer raw data}
\description{renal cell cancer raw data}
\usage{
data(qua)
data(spotanno)
data(hybanno)
}
\format{
\code{qua} is an \code{8704x4x175} array. Its first dimension are the spots
(probes) on the arrays, its second dimension are four different spot
statistics (\code{fg.green}, \code{bg.green}, \code{fg.red}, \code{bg.red}),
its third dimension are the 175 different arrays.

\code{spotanno} is a data frame with 8704 rows, one for each spot. Its
columns are:
\code{Block}, \code{Row}, \code{Column}, \code{AccNumber}, \code{Name},
\code{SrcPlt}, \code{SrcRow}, \code{SrcCol}, \code{vendor},  \code{ImageID},
\code{probe}. The column \code{spotanno$probe} is of class \code{numeric},
with values from 1 to 4224, referring to the rows of \code{\link{cloneanno}}.

\code{hybanno} is a data frame with 175 rows, one for each array. Its 
columns are:
\code{filename}, \code{patientid}, \code{slideid}. The column \code{patientid}
is of class \code{character}, refering to the patients in 
\code{pData(\link[kidpack:eset]{eset})}.
}

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

\seealso{\code{\link[kidpack:eset]{eset}}}
\keyword{datasets}
\examples{
data(qua)
qua[1,,1]
data(spotanno)
spotanno[1,]
data(hybanno)
hybanno[1,]

## see also the vignette
}

