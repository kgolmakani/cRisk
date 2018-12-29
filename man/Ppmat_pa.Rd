\name{Ppmat_pa}
\alias{Ppmat_pa}
\title{
The transition mtrix  for the last round when the total number of round attended is greater than 1.
}
\description{
Function for computing the transition mtrix  for the last round when the total number of round attended is greater than 1.
}
\usage{
Ppmat_pa (m1_pa, m2_pa)(h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{m1_pa}{
A vector of multinomial regression’s coffiecients considering the the last observations of subjects who attanded more than once and their current state is true negative. In this regression the result of screening is regressed on the censoring time.
}
\item{m2_pa}{
A vector of multinomial regression’s coffiecients considering the the last observations of subjects who attanded more than once and their current state is false positive. In this regression the result of screening is regressed on the censoring time.
}
\item{h}{
Total number of screening rounds attended. 
}
}

\value{
A transition matrix which it’s entries represent the transition probabilities where   1<t=h. The first row represent the transition probabilities when the initial state is true negative. The second row represent the transition probabilities when the initial state is false positive and the third row show the transition probabilities when the initial state is the competing risk. 
}
