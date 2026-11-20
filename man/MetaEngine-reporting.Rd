\name{MetaEngine-reporting}
\alias{MetaEngine-reporting}
\alias{pm_result_summary}
\alias{pm_print_summary}
\alias{pm_write_results_summary}
\title{MetaEngine reporting functions}
\description{Count differential features at specified nominal-p, FDR, and fold-change thresholds and write tabular and text summaries.}
\details{Significance thresholds use strict less-than comparisons. Fold-change thresholds use \code{2^abs(logFC)} with greater-than-or-equal comparisons. Up and down counts follow the sign of \code{logFC}, so their interpretation depends on the fitted contrast.}
\keyword{internal}
