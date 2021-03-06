% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_testing.R
\name{remove_recent}
\alias{remove_recent}
\title{Remove latest surveys by year.}
\usage{
remove_recent(df, threshold = 2015)
}
\arguments{
\item{df}{Input data frame.}

\item{threshold}{Cutoff year.}
}
\description{
Removes all surveys beyond the specified threshold year
}
\seealso{
Other testing functions: 
\code{\link{alt1_smpl2}()},
\code{\link{alt1_smpl3}()},
\code{\link{alt2_flat2}()},
\code{\link{alt2_flat3}()},
\code{\link{alt3_ltst3}()},
\code{\link{remove_latest}()},
\code{\link{test_abc2}()},
\code{\link{test_abc}()}
}
\concept{testing functions}
