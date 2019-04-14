\name{mucplot}
\alias{mucplot}
\title{Model Uncertainty Curve}
\description{
The function is to plot the model uncertainty curve based on the result of function mcb or mcb.compare.
}
\usage{
mucplot(obj)
}
\arguments{
\item{obj}{the result of function mcb or mcb.compare.}
}
\examples{
data(Diabetes) # load data
x <- Diabetes[,c('S1','S2','S3','S4','S5','S6')]
y <- Diabetes[,c('Y')]
x <- data.matrix(x)
y <- data.matrix(y)
result <- mcb(x=x, y=y)
# result <- mub.compare(x=x,y=y)
mucplot(result) # plot model uncertainty curve
}
