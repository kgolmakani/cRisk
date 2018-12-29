\name{Ppmat_nh}
\alias{Ppmat_nh}
\title{
The transition mtrix considering only the last round when subjects attend more than one round.
}
\description{
Function for computing the transition mtrix considering only the last round when subjects attend more than one round. The transition probabilities depend on screening round, censoring time and number of prior false positive at round g. 

}
\usage{
Ppmat_nh(m1,m2)(h,t,g,x)
}

\arguments{
\item{m1}{
A (2M)×5 matrix  where M is the maximum number of screening rounds and initial state is true negative. For i from 1 to M, m1_nh[2(i-1)+1,] is the coefficients of a multinomial regression when the state is false positive and m1_nh[2i+1,] is the coefficients of a multinomial regression when the state is competing event .In this multinomial regression the variable state is regressed on the number of prior false positive at round i. In this regression only last round of the subjects who attanded more than one round is included. 
}
\item{m2}{
A (2M)×5 matrix  where M is the maximum number of screening rounds and initial state is false positive. For i from 1 to M, m2_nh[2(i-1)+1,] is the coefficients of a multinomial regression when the state is  false positive and m2_nh[2i+1,] is the coefficients of a multinomial regression when the state is competing .In this multinomial regression the variable state is regressed on the number of prior false positive at round i. In this regression only last round of the subjects who attanded more than one round is included.
}
\item{h}{
vector of censoring time.
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
A transition matrix which it’s entries represent the transition probabilities where   1<t=h. The first row represent the transition probabilities when the initial state is true negative. The second row represent the transition probabilities when the initial state is false positive and the third row show the transition probabilities when the initial state is the competing risk. 
Transition probabilities depend on the number of prior false positives, censoring time, screening round and interaction between round and censoring time.
}
