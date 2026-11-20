\name{MetaEngine-preprocessing}
\alias{MetaEngine-preprocessing}
\alias{pm_load_data}
\alias{pm_filter_missing}
\alias{pm_impute_halfmin}
\alias{pm_transform}
\title{MetaEngine preprocessing functions}
\description{Read feature-intensity and sample-metadata tables, filter missing values, impute half-minimum intensities, and transform the resulting matrix.}
\details{Missingness filters count zero and non-finite intensities as missing. Imputation uses half the smallest positive finite value within each feature. Optional transformations include log2, quantile normalization, and feature-wise centering or scaling. These functions operate on processed feature tables.}
\keyword{internal}
