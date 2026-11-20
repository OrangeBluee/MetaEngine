# MetaEngine workflow module.

pm_run <- function(intensity,metadata,formula,contrast=NULL,
                   feature_id_col="mz_time",sample_id_col="sample_id",
                   annotation_cols=NULL,
                   max_missing=0.5,group_col=NULL,max_group_missing=NULL,
                   impute=TRUE,log2=TRUE,quantile=FALSE,
                   block_col=NULL,
                   p_thresholds=c(0.05,0.01,0.001),
                   fdr_thresholds=c(0.10,0.05,0.01),
                   fc_thresholds=c(1.5,2,3),
                   annotation_source=NULL,
                   annotation_file=NULL,
                   annotation_sheet=NULL,
                   annotation_library=NULL,
                   confirmed_annotations=NULL,
                   msms_evidence=NULL,
                   annotation_ppm=5,
                   annotation_rt_tol=0.30,
                   msms_score_threshold=0.70,
                   msdial_high_msms=0.80,
                   msdial_moderate_msms=0.60,
                   knowledge_base=NULL,
                   metabolite_database=NULL,
                   metabolite_database_type="auto",
                   metabolite_database_name=NULL,
                   metabolite_database_is_saliva=FALSE,
                   database_ppm=5,
                   ion_mode=NULL,
                   expand_adducts=FALSE,
                   empirical_compounds=TRUE,
                   empirical_cor_thresh=0.80,
                   empirical_rt_tol=0.15,
                   output_dir="MetaEngine_output") {

  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

  d <- pm_load_data(
    intensity, metadata,
    feature_id_col, sample_id_col,
    annotation_cols
  )

  n_samples_total <- nrow(d$metadata)

  if (!inherits(formula, "formula"))
    formula <- as.formula(formula)

  model_vars <- all.vars(formula)
  vars_to_check <- model_vars

  if (!is.null(block_col))
    vars_to_check <- unique(c(vars_to_check, block_col))

  missing_vars <- setdiff(vars_to_check, names(d$metadata))
  if (length(missing_vars)) {
    stop("Variables in model not found in metadata: ",
         paste(missing_vars, collapse=", "))
  }

  complete_idx <- complete.cases(
    d$metadata[, vars_to_check, drop=FALSE]
  )

  excluded_samples <- character(0)

  if (any(!complete_idx)) {
    excluded_samples <- as.character(
      d$metadata[[sample_id_col]][!complete_idx]
    )

    message("MetaEngine: excluding ",
            length(excluded_samples),
            " sample(s) with missing model covariates: ",
            paste(excluded_samples, collapse=", "))

    d$metadata <- d$metadata[complete_idx,,drop=FALSE]
    d$X <- d$X[,complete_idx,drop=FALSE]
  }

  d$metadata[] <- lapply(
    d$metadata,
    function(z) if (is.factor(z)) droplevels(z) else z
  )

  n_samples_used <- nrow(d$metadata)
  n_samples_excluded <- n_samples_total - n_samples_used

  grp <- if(!is.null(group_col)) d$metadata[[group_col]] else NULL

  X <- pm_filter_missing(
    d$X,
    max_missing,
    grp,
    max_group_missing
  )

  n_features_after_filtering <- nrow(X)

  if(impute)
    X <- pm_impute_halfmin(X)

  X <- pm_transform(
    X,
    log2=log2,
    quantile=quantile
  )

  # -----------------------------------------------------------------
  # Differential analysis
  # -----------------------------------------------------------------
  res <- pm_differential(
    X,
    d$metadata,
    formula,
    contrast,
    block_col
  )

  # Add ordinary fold change to the differential-results table.
  if ("logFC" %in% names(res)) {
    res$FC <- 2^res$logFC
    res$abs_FC <- 2^abs(res$logFC)
  }

  # -----------------------------------------------------------------
  # Optional annotation
  #
  # Direct MS-DIAL import
  #   annotation_source = "MSDIAL"
  #   annotation_file   = original MS-DIAL xlsx/txt/csv export
  #
  # Alternative input: local compound library with optional MS/MS evidence.
  # -----------------------------------------------------------------
  annotation_object <- NULL

  if (!is.null(annotation_source)) {

    source_key <- toupper(as.character(annotation_source)[1])

    if (source_key == "MSDIAL") {

      if (is.null(annotation_file))
        stop("annotation_file is required when annotation_source='MSDIAL'.")

      message("MetaEngine: importing MS-DIAL annotation evidence...")

      annotation_object <- pm_annotate_msdial(
        annotation_file=annotation_file,
        confirmed_annotations=confirmed_annotations,
        sheet=annotation_sheet,
        high_msms=msdial_high_msms,
        moderate_msms=msdial_moderate_msms
      )

    } else {
      stop(
        "Unsupported annotation_source: ", annotation_source,
        ". Currently supported: 'MSDIAL'."
      )
    }

  } else if (!is.null(annotation_library) || !is.null(confirmed_annotations)) {

    message("MetaEngine: running local-library metabolite annotation...")

    if (is.null(annotation_library)) {
      annotation_library_use <- data.frame(
        compound_id=character(),
        metabolite_name=character(),
        mz=numeric(),
        rt=numeric(),
        stringsAsFactors=FALSE
      )
    } else {
      annotation_library_use <- annotation_library
    }

    feature_info <- pm_parse_mz_time(rownames(X))

    annotation_object <- pm_annotate(
      features=feature_info,
      compound_library=annotation_library_use,
      confirmed_annotations=confirmed_annotations,
      msms_evidence=msms_evidence,
      ppm=annotation_ppm,
      rt_tol=annotation_rt_tol,
      msms_score_threshold=msms_score_threshold
    )
  }

  if (!is.null(annotation_object)) {
    # Restrict annotation to features retained after preprocessing/QC.
    keep_ann <- annotation_object$best$mz_time %in% rownames(X)
    annotation_object$best <- annotation_object$best[keep_ann,,drop=FALSE]

    if (!is.null(annotation_object$candidates) &&
        "mz_time" %in% names(annotation_object$candidates)) {
      annotation_object$candidates <- annotation_object$candidates[
        annotation_object$candidates$mz_time %in% rownames(X),,
        drop=FALSE
      ]
    }

    pm_write_annotation_outputs(
      annotation_object,
      output_dir
    )

    res <- pm_merge_annotation_results(
      res,
      annotation_object$best
    )
  }


  # -----------------------------------------------------------------
  # Optional saliva/oral/exposure knowledge-base integration
  # -----------------------------------------------------------------
  if (!is.null(annotation_object) && !is.null(knowledge_base)) {
    annotation_object$best <- pm_merge_knowledge_base(
      annotation_object$best,
      knowledge_base
    )

    annotation_object$best <- pm_annotation_priority(
      annotation_object$best
    )

    res <- pm_merge_annotation_results(
      res,
      annotation_object$best
    )

    write.csv(
      annotation_object$best,
      file.path(output_dir,"annotation_best_match_context_enriched.csv"),
      row.names=FALSE
    )
  }


  # -----------------------------------------------------------------
  # Independent metabolite-database annotation
  # -----------------------------------------------------------------
  database_annotation_object <- NULL

  if (!is.null(metabolite_database)) {

    if (is.null(ion_mode))
      stop("ion_mode is required when metabolite_database is supplied.")

    message(
      "MetaEngine: independently annotating all retained features against ",
      ifelse(
        is.null(metabolite_database_name),
        "metabolite database",
        metabolite_database_name
      ),
      "..."
    )

    feature_db <- pm_parse_mz_time(rownames(X))

    database_annotation_object <- pm_annotate_database(
      features=feature_db,
      database_file=metabolite_database,
      database_type=metabolite_database_type,
      database_name=metabolite_database_name,
      saliva_database=metabolite_database_is_saliva,
      ion_mode=tolower(ion_mode),
      ppm=database_ppm,
      msdial_best=if (is.null(annotation_object)) NULL else annotation_object$best,
      confirmed_annotations=confirmed_annotations
    )

    pm_write_database_annotation_outputs(
      database_annotation_object,
      output_dir
    )

    # Database-derived best candidate is the primary MetaEngine annotation.
    # MS-DIAL preliminary annotation remains in separate msdial_* columns.
    dbbest <- database_annotation_object$best

    db_cols <- intersect(
      c(
        "metabolite_name","compound_id","formula","inchikey","smiles",
        "kegg_id","pubchem_id","chebi_id","database",
        "saliva_database","adduct_hypothesis","neutral_mass_hypothesis",
        "database_ppm_error","database_annotation_score",
        "annotation_source","annotation_level","n_database_candidates",
        "ambiguous","msdial_preliminary_name","msdial_annotation_level",
        "msdial_annotation_score","msdial_name_agreement",
        "msdial_inchikey_agreement","msdial_formula_agreement",
        "msdial_msms_supported","confirmed","confirmation_source"
      ),
      names(dbbest)
    )

    ii <- match(res$mz_time, dbbest$mz_time)

    # Keep MS-DIAL-derived metabolite field distinguishable if it already exists.
    if ("metabolite_name" %in% names(res)) {
      names(res)[names(res)=="metabolite_name"] <- "msdial_metabolite_name"
    }
    if ("annotation_level" %in% names(res)) {
      names(res)[names(res)=="annotation_level"] <- "msdial_evidence_level"
    }
    if ("annotation_score" %in% names(res)) {
      names(res)[names(res)=="annotation_score"] <- "msdial_evidence_score"
    }

    for (nm in db_cols) {
      target_nm <- paste0("metaengine_", nm)
      res[[target_nm]] <- dbbest[[nm]][ii]
    }
  }

  # -----------------------------------------------------------------
  # Optional adduct expansion
  # -----------------------------------------------------------------
  adduct_expansion <- NULL

  if (isTRUE(expand_adducts)) {
    if (is.null(ion_mode))
      stop("ion_mode must be 'positive' or 'negative' when expand_adducts=TRUE.")

    fi <- pm_parse_mz_time(rownames(X))

    if (!is.null(annotation_object) &&
        "adduct" %in% names(annotation_object$best)) {
      fi$adduct <- annotation_object$best$adduct[
        match(fi$mz_time, annotation_object$best$mz_time)
      ]
    }

    adduct_expansion <- pm_expand_adducts(
      fi,
      mode=tolower(ion_mode),
      include_reported_adduct=TRUE
    )

    write.csv(
      adduct_expansion,
      file.path(output_dir,"adduct_expansion.csv"),
      row.names=FALSE
    )
  }

  # -----------------------------------------------------------------
  # Empirical compounds from RT proximity + correlation
  # -----------------------------------------------------------------
  empirical_map <- NULL

  if (isTRUE(empirical_compounds)) {
    empirical_map <- pm_empirical_compounds_mztime(
      X,
      cor_thresh=empirical_cor_thresh,
      rt_tol=empirical_rt_tol
    )

    write.csv(
      empirical_map,
      file.path(output_dir,"empirical_compounds.csv"),
      row.names=FALSE
    )

    if (!is.null(annotation_object)) {
      annotation_object$best$empirical_compound <-
        empirical_map$empirical_compound[
          match(annotation_object$best$mz_time, empirical_map$mz_time)
        ]

      res$empirical_compound <-
        empirical_map$empirical_compound[
          match(res$mz_time, empirical_map$mz_time)
        ]
    }
  }

  # Populate common annotation columns for downstream workflows.
  # Provenance-specific msdial_* and metaengine_* fields remain available.
  res <- pm_unify_annotation_results(res)

  # -----------------------------------------------------------------
  # Save core outputs
  # -----------------------------------------------------------------
  write.csv(
    res,
    file.path(output_dir,"differential_results.csv"),
    row.names=FALSE
  )

  write.csv(
    data.frame(mz_time=rownames(X), X, check.names=FALSE),
    file.path(output_dir,"processed_matrix.csv"),
    row.names=FALSE
  )

  write.csv(
    d$metadata,
    file.path(output_dir,"metadata_used.csv"),
    row.names=FALSE
  )

  if (length(excluded_samples)) {
    write.csv(
      data.frame(sample_id=excluded_samples),
      file.path(output_dir,"excluded_samples.csv"),
      row.names=FALSE
    )
  }

  sum_obj <- pm_write_results_summary(
    res,
    output_dir=output_dir,
    p_thresholds=p_thresholds,
    fdr_thresholds=fdr_thresholds,
    fc_thresholds=fc_thresholds,
    n_samples_total=n_samples_total,
    n_samples_used=n_samples_used,
    n_samples_excluded=n_samples_excluded,
    model=formula,
    contrast=contrast
  )

  analysis_overview <- data.frame(
    metric=c(
      "Original features",
      "Features after missing-value filtering",
      "Samples total",
      "Samples analyzed",
      "Samples excluded",
      "Features annotated"
    ),
    value=c(
      nrow(d$annotation),
      n_features_after_filtering,
      n_samples_total,
      n_samples_used,
      n_samples_excluded,
      if (is.null(annotation_object)) NA_integer_ else
        sum(annotation_object$best$annotation_level != "unannotated", na.rm=TRUE)
    )
  )

  write.csv(
    analysis_overview,
    file.path(output_dir,"analysis_overview.csv"),
    row.names=FALSE
  )

  pm_print_summary(
    sum_obj$summary,
    n_samples_total=n_samples_total,
    n_samples_used=n_samples_used,
    n_samples_excluded=n_samples_excluded,
    model=formula,
    contrast=contrast
  )

  if (!is.null(annotation_object)) {
    ann_n <- sum(
      annotation_object$best$annotation_level != "unannotated",
      na.rm=TRUE
    )
    message(
      "MetaEngine annotation: ",
      ann_n, " / ", nrow(annotation_object$best),
      " features received at least one annotation."
    )
  }

  saveRDS(
    list(
      X=X,
      metadata=d$metadata,
      annotation_input=d$annotation,
      annotations=annotation_object,
      database_annotation=database_annotation_object,
      database_annotation=database_annotation_object,
      adduct_expansion=adduct_expansion,
      empirical_compounds=empirical_map,
      stats=res,
      summary=sum_obj$summary,
      manuscript_summary=sum_obj$manuscript,
      excluded_samples=excluded_samples,
      settings=list(
        p_thresholds=p_thresholds,
        fdr_thresholds=fdr_thresholds,
        fc_thresholds=fc_thresholds,
        annotation_source=annotation_source,
        annotation_file=annotation_file,
        annotation_ppm=annotation_ppm,
        annotation_rt_tol=annotation_rt_tol,
        msms_score_threshold=msms_score_threshold,
        msdial_high_msms=msdial_high_msms,
        msdial_moderate_msms=msdial_moderate_msms,
        metabolite_database=metabolite_database,
        metabolite_database_type=metabolite_database_type,
        metabolite_database_name=metabolite_database_name,
        metabolite_database_is_saliva=metabolite_database_is_saliva,
        database_ppm=database_ppm,
        ion_mode=ion_mode,
        expand_adducts=expand_adducts,
        empirical_compounds=empirical_compounds,
        empirical_cor_thresh=empirical_cor_thresh,
        empirical_rt_tol=empirical_rt_tol
      )
    ),
    file.path(output_dir,"analysis_object.rds")
  )

  invisible(
    list(
      X=X,
      metadata=d$metadata,
      annotation=annotation_object,
      adduct_expansion=adduct_expansion,
      empirical_compounds=empirical_map,
      stats=res,
      summary=sum_obj$summary,
      manuscript_summary=sum_obj$manuscript,
      excluded_samples=excluded_samples
    )
  )
}

