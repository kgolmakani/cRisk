\name{disc.surv}
\alias{disc.surv}
\title{
Cumulative risk of multiple false positive screening test under discrete survival model. 
}
\description{
A function for calculating the cumulative risk of multiple false positive screening test under discrete survival model. 
}
\usage{
disc.surv (data,fpl.cens, M,l)
}

\arguments{
\item{data}{
A data frame that include S ( censoring time or total number of screening round attended).
}

\item{fpl.cens}{
A vector for the number of screening rounds of the l th false positive results that is either less than or greater than S or 0. 
}
\item{M}{
Maximum number of screening round.
}
\item{l}{
Number of false positive for computing cumulative risk.
}

}

\value{
A list which its first element is cumulative risk of  multiple false positive screening tests under discrete survival model and the second element is analytical standard error. 
}
