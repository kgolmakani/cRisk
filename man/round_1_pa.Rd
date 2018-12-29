\name{round_1_pa}
\alias{round_1_pa}
\title{
Probability of  each state in the first round. 
}
\description{
A function for computing the probability of  each state in the first round which depends on the censoring time.
}
\usage{
round_1_pa (m1_pa, m2_pa)(i,h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{m1_pa}{
A vector containing  the coefficients of multinomial regression cosidering subjects who attended only one round. In this regression the result of screening (states) is regressed on 1. 
}
\item{m2_pa}{
A vector containing the coefficients of logistic regression considering the first round of all subjects who attended more than one rounds. In this regression the result of screening (states) is regressed on total number of rounds attended (censoring time).
}
\item{i}{
State that can be 0 (true negative), 1 (false positive) or 2 ( competing event).
}
\item{h}{
total number of screening round attended.
}
}

\value{
A numeric value between 0 and 1 (probability) which represtns the probability of an specific state in the first round. 
}
