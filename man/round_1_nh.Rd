\name{round_1_nh}
\alias{round_1_nh}

\title{
Probability of each state (screening result) in the first round
}
\description{
Function for computing the probability of each state in the first round which depends on the censoring time and number of prior false positive results at round g. 
}
\usage{
round_1_nh(m1_nh, m2_nh)(i, h,g,x)
}

\arguments{
  \item{m1_nh}{
A (2M)×2 matrix  where M is the maximum number of screening rounds and initial state is true negative. For i from 1 to M, m1_nh[2(i-1)+1,] is the coefficients of a multinomial regression when the state is false positive and m1_nh[2i+1,] is the coefficients of a multinomial regression when the state is competing event.In this multinomial regression the variable state is regressed on the number of prior false positive at round i. In this regression subjects who attended only one round is included. 

}
\item{m2_nh}{
A (2M)×2 matrix  where M is the maximum number of screening rounds and initial state is false positive. For i from 1 to M, m2_nh[2(i-1)+1,] is the coefficients of a multinomial regression when the state is  false positive and m2_nh[2i+1,] is the coefficients of a multinomial regression when the state is competing event .In this multinomial regression the variable state is regressed on the number of prior false positive at round i. In this regression subjects who attended only one round is included. 
}
\item{h}{
Vector of total number of screening round attended.
}
\item{t}{
Vector of screening rounds.
}
\item{g}{
Numeric value representing a fixed screening round.
}
\item{x}{
A numeric value that represents the number of prior false positives at round g.
}
}

\value{
A numeric value between 0 and 1 (probability) which represtns the probability of an specific state in the first round.
}
