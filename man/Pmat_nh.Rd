\name{Pmat_nh}
\alias{Pmat_nh}

\title{
transition mtrix when the first and last round is not included.
}
\description{
function for computing the transition mtrix when the first and last round is not included and the transition probabilities depend on screening round, censoring time and number of prior false positive at round g. 
}
\usage{
Pmat_nh(m1,m2)(h,t,g,x)
}

\arguments{
  \item{m1}{
A M×5 matrix where M is the maximum number of screening round. For i from 1 to M, m1[i,] is the coefficients of a logistic regression for which the screening results (states) is regressed on number of prior false positives at round i, censoring time, round and interaction between round and censoring time. In this regression all observations of each subject (except their first and last observation) with false positive result are included.  
}
\item{m2}{
A M×5 matrix where M is the maximum number of screening round. For i from 1 to M, m1[i,] is the coefficients of a logistic regression for which the screening results (states) is regressed on number of prior false positives at round i, censoring time, round and interaction between round and censoring time. In this regression all observations of each subject (except their first and last observation) with false competing event are included.  
}
\item{h}{
vector of total number of screening rounds attended.
}
\item{t}{
vector of screening rounds.
}

\item{g}{
numeric value specifying a screening round.
}
\item{x}{
a numeric value that represents the number of prior false positives at round g.
}

}

\value{
A transition matrix which it’s entries represent the transition probabilities where   1<t<h. The first row represent the transition probabilities when the initial state is true negative. The second row represent the transition probabilities when the initial state is false positive and the third row show the transition probabilities when the initial state is the competing risk. Transition probabilities depend on the number of prior false positives, censoring time, screening round and interaction between round and censoring time.
}
