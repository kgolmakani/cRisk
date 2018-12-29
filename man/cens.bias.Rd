\name{cens.bias}
\alias{cens.bias}
\title{
Cumulative risk of multiple false positive screening tests based on censoring bias model.
}
\description{
Function for computing cumulative risk of multiple false positive screening tests based on censoring bias model.
}
\usage{
cens.bias (data, fpl, M, alpha, l, functions)
}

\arguments{
\item{data}{
A data frame that contains S (censoring time), StudyID_c (ID variable), round (screening round). 
}
\item{fpl}{
A vector of the round of the lth false positive for each observation of the data frame. 
}
\item{M}{
Maximum number of screening round.
}
\item{alpha}{
Censoring bias parameter.
}
\item{l}{
Number of false positive for computing cumulative risk.
}
\item{functions}{
A list contaning round_1(m1,m2) (function of censoring time and screening result which determine the probability of screening result at the first round), Ppmat (m1,m2) ( a transition matrix for the last round when the total number of round attended is greater than 1. each transition probability can be a function of censoring time, screening round,  interaction of both and patient characteristics), Pmat(m1,m2) (transition mtrix when the first and last round is not included. The transition probabilities depend on screening round, censoring time and interaction of both and patient charcateristics), round_1_nh(m1_nh, m2_nh) (probability of each state in the first round which depends on the censoring time and number of prior false positive results at round g. it can also be dependent on patient characteristics), Pmat_nh(m1,m2) (the transition mtrix when the first and last round is not included and the transition probabilities depend on screening round, censoring time and number of prior false positive at round g ans patient characteristics),  Ppmat_nh(m1,m2) (transition mtrix considering only the last round when subjects attend more than one round. The transition probabilities depend on screening round, censoring time and number of prior false positive at round g and patient characteristics).
}
}

\value{
A list that its first element is cumulative risk of l false positive screening round after M screening follow-up under censoring bias model. The second element is the analytical standard error. 
}
