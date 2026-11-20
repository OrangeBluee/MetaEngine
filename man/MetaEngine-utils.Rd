\name{MetaEngine-utils}
\alias{MetaEngine-utils}
\alias{pm_require}
\alias{pm_initialize}
\alias{pm_bh}
\alias{pm_read}
\title{MetaEngine input and dependency utilities}
\description{Read delimited tables, check package dependencies, initialize MetaEngine, and apply Benjamini-Hochberg adjustment to finite p-values.}
\details{Dependency installation is disabled by default and can be enabled with \code{install = TRUE}. Table input accepts a data frame or a local CSV or tab-delimited file. The p-value adjustment helper preserves non-finite entries as missing values.}
\keyword{internal}
