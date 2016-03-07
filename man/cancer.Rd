\name{cancer}
\alias{cancer}
\docType{data}
\title{
Cancer Data
}
\description{
Yearly cancer death rate (per 100,000 population) of Pennsylvania between 1930 and 2000.
}
\usage{data(cancer)}
\format{
  A time series with 71 observations on the following variable.
  \describe{
    \item{\code{cancer}}{a time series for yearly cancer death rate}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{

Wei, W. W.S. (2006). \emph{Time Series Analysis: Univariate and Multivariate Methods (2nd ed.)} Boston: Addison-Wesley.

}
\examples{
data(cancer)
plot(cancer)
out = fcast(cancer, f.number = 2)
out$mean
}
\keyword{datasets}
