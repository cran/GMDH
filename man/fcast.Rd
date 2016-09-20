\name{fcast}
\alias{fcast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Make Short Term Forecasting via GMDH-Type Neural Network Algorithms
}
\description{
\code{fcast} forecasts time series via GMDH-type neural network algorithms.
}

\usage{fcast(data, method = "GMDH", input = 4, layer = 3, f.number = 5, level = 95,
tf = "all", weight = 0.7,lambda = c(0,0.01,0.02,0.04,0.08,0.16,0.32,0.64,
1.28,2.56,5.12,10.24))}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{is a univariate time series of class ts}
  \item{method}{expects a character string to choose the desired method to forecast time series. To utilize GMDH-type neural network in forecasting, method is set to "GMDH". One should set method to "RGMDH" for forecasting via Revised GMDH-type neural network. Default is set to "GMDH"}
  \item{input}{is the number of inputs. Defaults input = 4} 
  \item{layer}{is the number of layers. Default is set to layer = 3} 
  \item{f.number}{is the number of observations to be forecasted. Defaults f.number = 5 }
  \item{level}{confidence level for prediction interval. Default is set to 95}
  \item{tf}{expects a character string to choose the desired transfer function to be used in forecasting. To use polynomial function, tf should be set to "polynomial". Similarly, tf should be set to "sigmoid", "RBF", "tangent" to utilize sigmoid function, radial basis function and tangent function, respectively. To use all functions simultaneously, default is set to "all"} 
  \item{weight}{is the percent of the data set to be utilized as learning set to estimate regularization parameter via validation. Default is set to weight = 0.70}
  \item{lambda}{is a vector which includes the sequence of feasible regularization parameters. Defaults lambda=c(0,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12,10.24)}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Returns a list containing following elements:
\item{method}{the forecasting method as a character string}
\item{mean}{point forecasts as a time series}
\item{lower}{lower limits for prediction interval}
\item{upper}{upper limits for prediction interval}
\item{level}{confidence level for prediction interval}
\item{x}{the original time series}
\item{fitted}{the fitted values}
\item{residuals}{the residuals of the model. The residuals are x minus the fitted values}

}

\references{

Dag, O., Yozgatligil, C. (2016). GMDH: An R Package for Short Term Forecasting via GMDH-Type Neural Network Algorithms. \emph{The R Journal}, \bold{8:1}, 379-386.


Ivakhnenko, A. G. (1966). Group Method of Data Handling - A Rival of the Method of Stochastic Approximation. \emph{Soviet Automatic Control}, \bold{13}, 43-71.

Kondo, T., Ueno, J. (2006). Revised GMDH-Type Neural Network Algorithm With A Feedback Loop Identifying Sigmoid Function Neural Network. \emph{International Journal of Innovative Computing, Information and Control}, \bold{2:5}, 985-996. 
}


\author{
Osman Dag, Ceylan Yozgatligil
}



\examples{
data = ts(rnorm(100, 10, 1))
out = fcast(data)
out

data = ts(rnorm(100, 10, 1))
out = fcast(data, input = 6, layer = 2, f.number = 1)
out$mean
out$fitted
out$residuals
plot(out$residuals)
hist(out$residuals)

}

\keyword{functions}
