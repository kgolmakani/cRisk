\name{round_1}
\alias{round_1}

\title{
Probability of  each state in the first round 
}
\description{
function for computing the probability of  each state in the first round which depends on the censoring time
}
\usage{
round_1(m1, m2)(i,h)
}

\arguments{
  \item{m1}{
vector of coefficients of multinomial regression for subjects who attended one screening round. In this regression variable states (screening results) is regressed on 1. 
}

\item{m2}{
vector of coeffiecient of logistic regression considering round 1 of all subjects who attanded more than one round. In this regression variable states (screening results) is regressed on censoring time. 
}

\item{i}{
state which can be 0 (true negative), 1 (false positive) or 2 ( competing event).
}

\item{h}{
total number of screening round attended.
}
}

\value{
A numeric value between 0 and 1 which represents the probability of a specific state in the first round.
}
