\name{fcast}
\alias{fcast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Predict and Forecast Time Series via GMDH-Type Neural Network Algorithms
}
\description{
\code{fcast} predicts and forecasts time series via GMDH-type neural network algorithms.
}

\usage{fcast(data, method = "GMDH", input = 4, layer = 3, f.number = 10, tf = "all", 
plotit = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{is an univariate time series}
  \item{method}{expects a character string to choose the desired method to forecast time series. To utilize GMDH-type neural network in forecasting, method is set to "GMDH". One should set method to "RGMDH" for forecasting via Revised GMDH-type neural network. Default is set to "GMDH"}
  \item{input}{is the number of inputs. Defaults input = 4} 
  \item{layer}{is the number of layers. Default is set to layer = 3} 
  \item{f.number}{is the number of observations to be forecasted. Defaults f.number = 10 }
  \item{tf}{expects a character string to choose the desired transfer function to be used in forecasting. To use polynomial function, tf should be set to "polynomial". Similarly, tf should be set to "sigmoid", "RBF", "tangent" to utilize sigmoid function, radial basis function and tangent function, respectively. To use all functions simultaneously, default is set to "all"} 
  \item{plotit}{is logical which controls whether historical data with forecasts should be plotted. Defaults plotit = TRUE}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Returns a list containing following elements:
\item{fitted}{fitted values}
\item{MSE}{MSE of prediction}
\item{forecasts}{forecasts }

}

\references{

Dag, O., Yozgatligil, C. (2015). GMDH: An R Package for Predicting and Forecasting Time Series via GMDH-Type Neural Network Algorithms. To be submitted.

Ivakhnenko, A. G. (1966). Group Method of Data Handling - A Rival of the Method of Stochastic Approximation. \emph{Soviet Automatic Control}, \bold{13}, 43-71.

Kondo, T., Ueno, J. (2006). Revised GMDH-Type Neural Network Algorithm With A Feedback Loop Identifying Sigmoid Function Neural Network. \emph{International Journal of Innovative Computing, Information and Control}, \bold{2:5}, 985-996. 
}


\author{
Osman Dag, Ceylan Yozgatligil
}

\note{
This is the version 1.0 of this user documentation file.
}


\examples{
data = rnorm(100, 10, 1)
out = fcast(data)
out

data = rnorm(100, 10, 1)
out = fcast(data, input = 6, layer = 2, f.number = 5)
out$forecasts
out$fitted
out$MSE

}

\keyword{functions}
