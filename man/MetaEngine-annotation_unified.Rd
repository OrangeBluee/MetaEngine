\name{MetaEngine-annotation_unified}
\alias{MetaEngine-annotation_unified}
\alias{pm_unify_annotation_results}
\title{Unified annotation result fields}
\description{Populate common metabolite-name, source, level, and score fields while retaining source-specific annotation columns.}
\details{A named MetaEngine database candidate supplies the common fields when available; otherwise a named MS-DIAL result is used. This source precedence does not compare identification confidence across engines. Scores from different engines are heuristic and are not calibrated probabilities.}
\keyword{internal}
