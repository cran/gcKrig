\name{plot.mlegc}
\alias{plot.mlegc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Geostatistical Data and Fitted Mean
}

\usage{
\method{plot}{mlegc}(x, xlab = "xloc", ylab = "yloc", xlim = NULL, ylim = NULL,
     displaymain = FALSE, pch = 20, textcex = 0.8, plotcex = 1, angle = 60,
     col = 4, col.regions = gray(90:0/100),...)
          }

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{mlegc} inherited from function \code{\link{mlegc}}.
}

  \item{xlab, ylab}{a title for the x and y axis.}

  \item{xlim, ylim}{numeric vectors of length 2, giving the x and y coordinates ranges.
  if they equal to \code{NULL} then they will be adjusted from the data.}

  \item{displaymain}{if \code{displaymain = FALSE}, no overall title will be added to each plot.
  If \code{displaymain = TRUE}, some titles will be added accordingly.}

  \item{pch}{plotting character, i.e., the symbol to use in the 3-D scatter plot.}

  \item{textcex}{a numerical value giving the amount by which plotting text should be magnified relative to the default.}

    \item{plotcex}{a numerical value giving the amount by which plotting symbols should be magnified relative to the default.}

  \item{angle}{angle between x and y axis.}

  \item{col}{color of the text.}

  \item{col.regions}{color vector to be used reflecting magnitude of the dataset at different locations.
  The general idea is that this should be a color vector of gradually varying color.}

  \item{...}{further arguments passed to plot and panel settings.}



}
\description{
Four plots will be generated. A level plot with the number of counts at given locations; a
level plot with fitted values; a 3-D scatter plot with the original counts; a
3-D scatter plot with fitted values.
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.simgc}},
\code{\link{plot.predgc}}
}


\keyword{Data Visualization}
\keyword{Estimation}
