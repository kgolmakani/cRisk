\name{Ppmat}
\alias{Ppmat}

\title{
Transition mtrix  for the last round when the total number of round attended is greater than 1.
}
\description{
function for computing the transition mtrix  for the last round when the total number of round attended is greater than 1.
}
\usage{
Ppmat(m1 , m2)(h,t)
}

\arguments{
  \item{m1}{
vector of coefficients of  multinomial regression when states (screening round results) are regressed on  screening round, censoring time and interaction of  both.  In this regression only last screening round  is considered when total number of rounds attended is greater than 1 and true negative is the current state. 
}
\item{m2}{
vector of coefficients of  multinomial regression when states (screening round results) are regressed on  screening round, censoring time and interaction of  both.  In this regression only last screening round  is considered when total number of rounds attended is greater than 1 and false positive is the current state.
}
\item{h}{
total number of screening round attended.
}
\item{t}{
screening round.
}
}

\value{
A transition matrix which it’s entries represent the transition probabilities where 1<t=h. The first row represent the transition probabilities when the initial state is true negative. The second row represent the transition probabilities when the initial state is false positive and the third row show the transition probabilities when the initial state is the competing risk. 
}
