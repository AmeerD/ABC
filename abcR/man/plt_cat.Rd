% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_plotting.R
\name{plt_cat}
\alias{plt_cat}
\title{Plot ABC Results by Category}
\usage{
plt_cat(
  mdl_proj,
  mdl_obs,
  level = "prim",
  sexes = "total",
  target_yr = 2020,
  countries = sample(unique(mdl_obs$country), 3)
)
}
\arguments{
\item{mdl_proj}{Projected values data frame.}

\item{mdl_obs}{Observed values data frame.}

\item{level}{Level to plot.}

\item{sexes}{Sex to plot.}

\item{target_yr}{Final year to plot.}

\item{countries}{Vector of countries to plot.}
}
\description{
\code{plt_cat} creates a faceted plot displaying the observed and
estimated true completion rates for a given level/sex combination.
Each facet presents a different country. The choice of countries can
be user provided, or random based on the observations provided.
}
\seealso{
Other plotting functions: 
\code{\link{mk_bias_plot}()},
\code{\link{mk_input_plot}()},
\code{\link{mk_input_survey_plot}()},
\code{\link{mk_par_plot}()},
\code{\link{mk_reg_plot}()},
\code{\link{mk_repl_plot}()},
\code{\link{mk_rhats_plot}()},
\code{\link{mk_trace_plot}()},
\code{\link{plt_country}()},
\code{\link{se_plotter}()}
}
\concept{plotting functions}
