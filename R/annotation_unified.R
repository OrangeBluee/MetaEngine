# Unified public annotation fields. Source-specific evidence remains untouched.
pm_unify_annotation_results <- function(results) {
  x <- results
  n <- nrow(x)
  get <- function(nm) if (nm %in% names(x)) as.character(x[[nm]]) else rep(NA_character_, n)
  valid <- function(z) !is.na(z) & nzchar(trimws(z)) &
    !tolower(trimws(z)) %in% c("unknown", "unidentified", "no match")
  db <- get("metaengine_metabolite_name")
  md <- if ("msdial_metabolite_name" %in% names(x)) get("msdial_metabolite_name") else get("metabolite_name")
  use_db <- valid(db)
  use_md <- !use_db & valid(md)
  x$metabolite_name <- ifelse(use_db, db, ifelse(use_md, md, NA_character_))
  x$annotation_source <- ifelse(use_db, "MetaEngine_database", ifelse(use_md, "MSDIAL", NA_character_))
  db_level <- get("metaengine_annotation_level")
  md_level <- if ("msdial_evidence_level" %in% names(x)) get("msdial_evidence_level") else get("annotation_level")
  x$annotation_level <- ifelse(use_db, db_level, ifelse(use_md, md_level, "unannotated"))
  db_score <- if ("metaengine_database_annotation_score" %in% names(x)) as.numeric(x$metaengine_database_annotation_score) else rep(NA_real_,n)
  md_score <- if ("msdial_evidence_score" %in% names(x)) as.numeric(x$msdial_evidence_score) else if ("annotation_score" %in% names(x)) as.numeric(x$annotation_score) else rep(NA_real_,n)
  x$annotation_score <- ifelse(use_db, db_score, ifelse(use_md, md_score, NA_real_))
  # Scores from distinct engines are heuristic and not calibrated probabilities.
  x
}
