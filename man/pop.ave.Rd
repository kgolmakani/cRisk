\name{pop.avg}
\alias{pop.avg}

\title{
Cumulative risk of multiple false positive screening test under population average model. 
}
\description{
A function for calculating the cumulative risk of multiple false positive screening test under population average model. 
}
\usage{
pop.avg (data, fpl.cens, functions1, functions2,l,M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A data frame containing S (censoring time), StudyID_c (ID), comp (the round when competing event occurs), delta (censoring status).
}
\item{fpl.cens}{
A vector for the number of screening rounds of the l th false positive results that is either less than or greater than S or 0. 
}
\item{functions1}{
A list contaning round_1_pa(m1,m2) (function of censoring time which determine the probability of screening result at the first round), Ppmat_pa (m1,m2) ( a transition matrix for the last round when the total number of round attended is greater than 1. each transition probability can be a function of censoring time and patient characteristics), Pmat_pa(m1,m2) (transition mtrix when the first and last round is not included. The transition probabilities depend on censoring time and interaction of both and patient charcateristics).
}
\item{functions2}{
A list contaning round_1(m1,m2) (function of screening round and censoring time which determine the probability of screening result at the first round), Ppmat (m1,m2) ( a transition matrix for the last round when the total number of round attended is greater than 1. each transition probability can be a function of censoring time, screening round,  interaction of both and patient characteristics), Pmat(m1,m2) (transition mtrix when the first and last round is not included. The transition probabilities depend on screening round, censoring time and interaction of both and patient charcateristics). 
}
\item{l}{
Number of false positive for computing cumulative risk.
}
\item{M}{
Maximum number of screening roud.
}
}
