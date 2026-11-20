# MetaEngine knowledge-base module.

pm_read_knowledge_base <- function(x) {
  if (is.null(x)) return(NULL)
  kb <- pm_read(x)

  if (!"metabolite_name" %in% names(kb))
    stop("Knowledge-base file requires metabolite_name.")

  kb$metabolite_name <- as.character(kb$metabolite_name)

  # Normalize supported Boolean columns when present.
  bool_cols <- intersect(
    c(
      "saliva_evidence",
      "oral_microbiome_related",
      "periodontal_related",
      "tobacco_related",
      "vaping_related",
      "oxidative_stress_related",
      "inflammation_related",
      "host_derived",
      "microbial_derived"
    ),
    names(kb)
  )

  for (nm in bool_cols)
    kb[[nm]] <- pm_flag(kb[[nm]])

  kb
}

pm_merge_knowledge_base <- function(annotation_best, knowledge_base) {
  if (is.null(knowledge_base)) return(annotation_best)

  an <- pm_read(annotation_best)
  kb <- pm_read_knowledge_base(knowledge_base)

  if (!"metabolite_name" %in% names(an))
    return(an)

  key_an <- tolower(trimws(as.character(an$metabolite_name)))
  key_kb <- tolower(trimws(as.character(kb$metabolite_name)))
  idx <- match(key_an, key_kb)

  extra <- setdiff(names(kb), "metabolite_name")

  for (nm in extra)
    an[[nm]] <- kb[[nm]][idx]

  an
}

pm_annotation_priority <- function(annotation_best) {
  x <- annotation_best

  evidence <- if ("annotation_score" %in% names(x)) {
    ifelse(is.finite(x$annotation_score), x$annotation_score, 0)
  } else rep(0, nrow(x))

  saliva <- if ("saliva_evidence" %in% names(x))
    ifelse(x$saliva_evidence %in% TRUE, 1, 0) else 0

  oral <- if ("oral_microbiome_related" %in% names(x))
    ifelse(x$oral_microbiome_related %in% TRUE, 1, 0) else 0

  exposure <- rep(0, nrow(x))
  for (nm in intersect(
    c("tobacco_related","vaping_related","periodontal_related"),
    names(x)
  )) {
    exposure <- pmax(
      exposure,
      ifelse(x[[nm]] %in% TRUE, 1, 0)
    )
  }

  # Calculate context and priority scores separately from identification evidence.
  x$context_score <- (
    0.50 * saliva +
    0.25 * oral +
    0.25 * exposure
  )

  x$priority_score <- (
    0.80 * evidence +
    0.20 * x$context_score
  )

  x
}

pm_write_knowledge_base_template <- function(output_dir=".") {
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

  kb <- data.frame(
    metabolite_name=c("Cotinine","Glutamate"),
    hmdb_id=c("HMDB0001046","HMDB0000148"),
    kegg_id=c("C00745","C00025"),
    saliva_evidence=c(TRUE,TRUE),
    oral_microbiome_related=c(FALSE,TRUE),
    periodontal_related=c(FALSE,TRUE),
    tobacco_related=c(TRUE,FALSE),
    vaping_related=c(TRUE,FALSE),
    oxidative_stress_related=c(FALSE,TRUE),
    inflammation_related=c(FALSE,TRUE),
    host_derived=c(TRUE,TRUE),
    microbial_derived=c(FALSE,TRUE),
    evidence_note=c(
      "Example exposure biomarker",
      "Example host/microbial metabolite"
    ),
    stringsAsFactors=FALSE
  )

  write.csv(
    kb,
    file.path(output_dir,"MetaEngine_knowledge_base_TEMPLATE.csv"),
    row.names=FALSE
  )

  invisible(kb)
}



