\name{Pmat_pa}
\alias{Pmat_pa}
\title{
Transition mtrix when the first and last round is not included.
}
\description{
Function for computing the transition mtrix when the first and last round is not included. The transition probabilities depend only on censoring time.      
}
\usage{
Pmat_pa (m1_pa , m2_pa) (h)
}

\arguments{
\item{m1_pa}{
A vector of  logistic regression coefficients considering only observation with current state true negative. In this regression results of screening is regressed on the censoring time.
}
\item{m2_pa}{
A vector of  logistic regression coefficients considering only observation with current state false positive. In this regression results of screening is regressed on the censoring time. 
}
\item{h}{
Total number of screening rounds attended.
}
}

\value{
A transition matrix which it’s entries represent the transition probabilities where 1<t<h. The first row represent the transition probabilities when the initial state is true negative. The second row represent the transition probabilities when the initial state is false positive and the third row show the transition probabilities when the initial state is the competing risk. 
}
