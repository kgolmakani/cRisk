\name{Pmat}
\alias{Pmat}

\title{
Transition mtrix when the first and last round is not included.
}
\description{
function for computing the transition mtrix when the first and last round is not included. The transition probabilities depend on screening round, censoring time and interaction of both.       
}
\usage{
Pmat(m1 , m2)(h,t)
}

\arguments{
  \item{m1}{
vector of coefficients of logistic regression considering all observations excluding first and last round with true negative as the current state . In this regression variable states (screening results) is regressed on screening round, censoring time and interaction of both. 
}
\item{m2}{
vector of coefficients of logistic regression considering all observations excluding first and last round with false positive as the current state . In this regression variable states (screening results) is regressed on screening round, censoring time and interaction of both.
}
\item{h}{
total number of screening round attended.
}
\item{t}{
screening round.
}
}

\value{
A transition matrix which it’s entries represent the transition probabilities where 1<t<h. The first row represent the transition probabilities when the initial state is true negative. The second row represent the transition probabilities when the initial state is false positive and the third row show the transition probabilities when the initial state is the competing risk. 
}
