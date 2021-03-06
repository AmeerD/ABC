% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_testing.R
\name{test_abc2}
\alias{test_abc2}
\title{Survey level ABC test.}
\usage{
test_abc2(input, tests, mcmc)
}
\arguments{
\item{input}{Input set data frame.}

\item{tests}{Test set data frame.}

\item{mcmc}{Raw ABC MCMC output.}
}
\value{
Returns a summary table with the model performance.
}
\description{
\code{test_abc2} uses the results of the ABC model executed using the input
set to estimate the completion rate indicator of the observation year of
each survey in the test set. The differences are summarised by the mean
squared error, mean absolute deviation, and median absolute deviation.
}
\seealso{
Other testing functions: 
\code{\link{alt1_smpl2}()},
\code{\link{alt1_smpl3}()},
\code{\link{alt2_flat2}()},
\code{\link{alt2_flat3}()},
\code{\link{alt3_ltst3}()},
\code{\link{remove_latest}()},
\code{\link{remove_recent}()},
\code{\link{test_abc}()}
}
\concept{testing functions}
