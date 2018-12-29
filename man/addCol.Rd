\name{addCol}
\alias{addCol}

\title{
Previous value of a longitudinal variable
}
\description{
add the variable of previous value of the variable state to the data frame. The data frame should include studyID_c ( id variable), round (screening round), result (current state or screening result).

}
\usage{
addCol(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
A data frame that includes studyID_c ( id variable), round (screening round), result (current state or screening result).
}
}

\value{
A data frame which also includes a variable for previous state or result.
}
