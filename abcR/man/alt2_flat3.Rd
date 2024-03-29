% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_deprecated_testing.R
\name{alt2_flat3}
\alias{alt2_flat3}
\title{Survey level flat model test.}
\usage{
alt2_flat3(input, tests)
}
\arguments{
\item{input}{Input set data frame.}

\item{tests}{Test set data frame.}
}
\value{
Returns a summary table with the model performance.
}
\description{
\code{alt2_flat3} uses the results of an intercept only model in the probit
space estimated using the input set to estimate the completion rate indicator
of the observation year of each survey in the test set. The differences are
summarised by the mean squared error, mean absolute deviation, and median
absolute deviation.
}
\seealso{
Other testing functions: 
\code{\link{alt1_smpl2}()},
\code{\link{alt1_smpl3}()},
\code{\link{alt2_flat2}()},
\code{\link{alt3_ltst3}()},
\code{\link{remove_latest}()},
\code{\link{remove_recent}()},
\code{\link{test_abc2}()},
\code{\link{test_abc}()}
}
\concept{testing functions}
