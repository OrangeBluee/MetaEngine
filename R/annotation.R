# MetaEngine annotation module.
# Database candidate matching, MS-DIAL evidence import, and annotation reporting.
# Candidate tables retain alternative identities and their supporting evidence.

pm_write_annotation_templates <- function(output_dir=".") {

  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

  compound_library <- data.frame(
    compound_id=c("HMDB0000122","HMDB0000056"),
    metabolite_name=c("Example metabolite A","Example metabolite B"),
    mz=c(150.0584, 180.0634),
    rt=c(2.315, 5.740),
    formula=c("C6H6N4O2","C6H12O6"),
    adduct=c("[M+H]+","[M-H]-"),
    database=c("HMDB","HMDB"),
    stringsAsFactors=FALSE
  )

  confirmed <- data.frame(
    mz_time=c("150.05840_2.315"),
    compound_id=c("HMDB0000122"),
    metabolite_name=c("Example metabolite A"),
    source=c("Authentic standard/internal standard"),
    stringsAsFactors=FALSE
  )

  msms <- data.frame(
    mz_time=c("180.06340_5.740"),
    compound_id=c("HMDB0000056"),
    spectral_score=c(0.91),
    source=c("Local MS/MS library"),
    stringsAsFactors=FALSE
  )

  write.csv(
    compound_library,
    file.path(output_dir,"annotation_compound_library_TEMPLATE.csv"),
    row.names=FALSE
  )

  write.csv(
    confirmed,
    file.path(output_dir,"annotation_confirmed_TEMPLATE.csv"),
    row.names=FALSE
  )

  write.csv(
    msms,
    file.path(output_dir,"annotation_msms_TEMPLATE.csv"),
    row.names=FALSE
  )

  invisible(
    list(
      compound_library=compound_library,
      confirmed=confirmed,
      msms=msms
    )
  )
}

pm_write_msdial_annotation_templates <- function(output_dir=".") {

  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

  confirmed <- data.frame(
    mz_time=c("151.06230_3.520"),
    metabolite_name=c("Example confirmed metabolite"),
    compound_id=c("optional_database_id"),
    source=c("Authentic standard/internal standard"),
    confirmation_method=c("m/z + RT + MS/MS"),
    notes=c("Replace example with validated metabolites"),
    stringsAsFactors=FALSE
  )

  write.csv(
    confirmed,
    file.path(output_dir,"confirmed_annotations_TEMPLATE.csv"),
    row.names=FALSE
  )

  invisible(confirmed)
}

# -----------------------------------------------------------------------------
# DATA
# -----------------------------------------------------------------------------

pm_ppm <- function(obs, theo) 1e6 * (obs-theo)/theo

pm_parse_mz_time <- function(mz_time) {
  x <- as.character(mz_time)

  # Parse m/z and retention time at the final underscore, e.g. 301.21668_8.417.
  pos <- regexpr("_[^_]*$", x)

  mz <- rep(NA_real_, length(x))
  rt <- rep(NA_real_, length(x))

  good <- pos > 1

  if (any(good)) {
    mz[good] <- suppressWarnings(
      as.numeric(substr(x[good], 1, pos[good]-1))
    )
    rt[good] <- suppressWarnings(
      as.numeric(substr(x[good], pos[good]+1, nchar(x[good])))
    )
  }

  data.frame(
    mz_time = x,
    mz = mz,
    rt = rt,
    stringsAsFactors = FALSE
  )
}

pm_annotation_features <- function(x,
                                   mz_time_col = "mz_time",
                                   mz_col = "mz",
                                   rt_col = "rt") {

  if (is.matrix(x)) {
    ids <- rownames(x)
    if (is.null(ids))
      stop("Intensity matrix must have mz_time row names for annotation.")

    return(pm_parse_mz_time(ids))
  }

  z <- pm_read(x)

  if (!mz_time_col %in% names(z))
    stop("Annotation feature table requires column: ", mz_time_col)

  out <- data.frame(
    mz_time = as.character(z[[mz_time_col]]),
    stringsAsFactors = FALSE
  )

  if (mz_col %in% names(z)) {
    out$mz <- suppressWarnings(as.numeric(z[[mz_col]]))
  } else {
    out$mz <- pm_parse_mz_time(out$mz_time)$mz
  }

  if (rt_col %in% names(z)) {
    out$rt <- suppressWarnings(as.numeric(z[[rt_col]]))
  } else {
    out$rt <- pm_parse_mz_time(out$mz_time)$rt
  }

  out
}

pm_standardize_library <- function(compounds,
                                   compound_id_col = "compound_id",
                                   name_col = "metabolite_name",
                                   mz_col = "mz",
                                   rt_col = "rt") {

  cdb <- pm_read(compounds)

  required <- c(compound_id_col, mz_col)
  miss <- setdiff(required, names(cdb))

  if (length(miss)) {
    stop(
      "Annotation library is missing required column(s): ",
      paste(miss, collapse=", ")
    )
  }

  out <- cdb

  names(out)[names(out) == compound_id_col] <- "compound_id"
  names(out)[names(out) == mz_col] <- "library_mz"

  if (name_col %in% names(out))
    names(out)[names(out) == name_col] <- "metabolite_name"

  if (!"metabolite_name" %in% names(out))
    out$metabolite_name <- as.character(out$compound_id)

  if (rt_col %in% names(out))
    names(out)[names(out) == rt_col] <- "library_rt"
  else
    out$library_rt <- NA_real_

  out$compound_id <- as.character(out$compound_id)
  out$metabolite_name <- as.character(out$metabolite_name)
  out$library_mz <- suppressWarnings(as.numeric(out$library_mz))
  out$library_rt <- suppressWarnings(as.numeric(out$library_rt))

  out
}

pm_match_mass <- function(features,
                          compounds,
                          ppm = 5,
                          rt_tol = NULL,
                          mz_time_col = "mz_time",
                          feature_mz = "mz",
                          feature_rt = "rt",
                          compound_id_col = "compound_id",
                          compound_name_col = "metabolite_name",
                          compound_mz = "mz",
                          compound_rt = "rt") {

  f <- pm_annotation_features(
    features,
    mz_time_col = mz_time_col,
    mz_col = feature_mz,
    rt_col = feature_rt
  )

  cdb <- pm_standardize_library(
    compounds,
    compound_id_col = compound_id_col,
    name_col = compound_name_col,
    mz_col = compound_mz,
    rt_col = compound_rt
  )

  out <- vector("list", nrow(f))

  for (i in seq_len(nrow(f))) {

    obs_mz <- f$mz[i]

    if (!is.finite(obs_mz))
      next

    dppm <- abs(pm_ppm(obs_mz, cdb$library_mz))

    hit <- which(
      is.finite(dppm) &
      dppm <= ppm
    )

    if (!length(hit))
      next

    z <- cdb[hit,,drop=FALSE]

    z$mz_time <- f$mz_time[i]
    z$observed_mz <- obs_mz
    z$observed_rt <- f$rt[i]
    z$ppm_error <- pm_ppm(obs_mz, z$library_mz)
    z$rt_error <- abs(z$observed_rt - z$library_rt)

    # Retain mass candidates and record RT agreement as a separate evidence field.
    z$rt_matched <- FALSE

    if (!is.null(rt_tol)) {
      z$rt_matched <- (
        is.finite(z$observed_rt) &
        is.finite(z$library_rt) &
        z$rt_error <= rt_tol
      )
    }

    out[[i]] <- z
  }

  ans <- do.call(rbind, out)

  if (is.null(ans)) {
    return(
      data.frame(
        mz_time=character(),
        compound_id=character(),
        metabolite_name=character(),
        observed_mz=numeric(),
        observed_rt=numeric(),
        library_mz=numeric(),
        library_rt=numeric(),
        ppm_error=numeric(),
        rt_error=numeric(),
        rt_matched=logical(),
        stringsAsFactors=FALSE
      )
    )
  }

  rownames(ans) <- NULL
  ans
}

pm_read_confirmed_annotations <- function(x) {

  if (is.null(x))
    return(NULL)

  z <- pm_read(x)

  required <- c("mz_time", "metabolite_name")
  miss <- setdiff(required, names(z))

  if (length(miss)) {
    stop(
      "Confirmed-annotation file requires: ",
      paste(required, collapse=", ")
    )
  }

  if (!"compound_id" %in% names(z))
    z$compound_id <- z$metabolite_name

  z$mz_time <- as.character(z$mz_time)
  z$compound_id <- as.character(z$compound_id)
  z$metabolite_name <- as.character(z$metabolite_name)

  z
}

pm_read_msms_evidence <- function(x) {

  if (is.null(x))
    return(NULL)

  z <- pm_read(x)

  required <- c("mz_time", "compound_id", "spectral_score")
  miss <- setdiff(required, names(z))

  if (length(miss)) {
    stop(
      "MS/MS evidence file requires columns: ",
      paste(required, collapse=", ")
    )
  }

  z$mz_time <- as.character(z$mz_time)
  z$compound_id <- as.character(z$compound_id)
  z$spectral_score <- suppressWarnings(as.numeric(z$spectral_score))

  z
}

pm_score_annotations <- function(candidates,
                                 confirmed_annotations = NULL,
                                 msms_evidence = NULL,
                                 ppm = 5,
                                 rt_tol = 0.30,
                                 msms_score_threshold = 0.70) {

  x <- candidates

  if (is.null(x) || !nrow(x))
    return(x)

  # ---------------------------------------------------------------
  # Join optional MS/MS evidence
  # ---------------------------------------------------------------
  x$spectral_score <- NA_real_
  x$msms_source <- NA_character_

  msms <- pm_read_msms_evidence(msms_evidence)

  if (!is.null(msms)) {

    key_x <- paste(x$mz_time, x$compound_id, sep="||")
    key_m <- paste(msms$mz_time, msms$compound_id, sep="||")

    idx <- match(key_x, key_m)

    x$spectral_score <- msms$spectral_score[idx]

    if ("source" %in% names(msms))
      x$msms_source <- as.character(msms$source[idx])
  }

  # ---------------------------------------------------------------
  # Quantitative evidence scores
  # ---------------------------------------------------------------
  x$mass_score <- pmax(
    0,
    1 - abs(x$ppm_error) / ppm
  )

  x$rt_score <- 0

  have_rt <- (
    is.finite(x$observed_rt) &
    is.finite(x$library_rt)
  )

  if (is.finite(rt_tol) && rt_tol > 0) {
    x$rt_score[have_rt] <- pmax(
      0,
      1 - x$rt_error[have_rt] / rt_tol
    )
  }

  x$msms_score <- ifelse(
    is.finite(x$spectral_score),
    pmin(pmax(x$spectral_score, 0), 1),
    0
  )

  # Combine mass, RT, and MS/MS scores; missing RT or MS/MS contributes zero.
  x$annotation_score <- (
    0.50 * x$mass_score +
    0.20 * x$rt_score +
    0.30 * x$msms_score
  )

  # ---------------------------------------------------------------
  # Annotation evidence levels
  # ---------------------------------------------------------------
  x$annotation_level <- "4_mass_only"

  x$annotation_level[
    isTRUE(rt_tol > 0) &
    x$rt_matched
  ] <- "3_mass_rt"

  x$annotation_level[
    is.finite(x$spectral_score) &
    x$spectral_score >= msms_score_threshold
  ] <- "2_msms_library"

  x$confirmed <- FALSE
  x$confirmation_source <- NA_character_

  conf <- pm_read_confirmed_annotations(confirmed_annotations)

  if (!is.null(conf)) {

    # Exact mz_time confirmation overrides computational candidates.
    conf_key <- conf$mz_time

    for (i in seq_len(nrow(conf))) {

      hit <- which(x$mz_time == conf$mz_time[i])

      if (!length(hit)) {

        # Add confirmed annotation even when absent from reference library.
        newrow <- x[1,,drop=FALSE]
        newrow[,] <- NA

        newrow$mz_time <- conf$mz_time[i]
        newrow$compound_id <- conf$compound_id[i]
        newrow$metabolite_name <- conf$metabolite_name[i]

        parsed <- pm_parse_mz_time(conf$mz_time[i])
        newrow$observed_mz <- parsed$mz
        newrow$observed_rt <- parsed$rt

        newrow$annotation_score <- 1
        newrow$annotation_level <- "1_confirmed_standard"
        newrow$confirmed <- TRUE

        if ("source" %in% names(conf))
          newrow$confirmation_source <- as.character(conf$source[i])
        else
          newrow$confirmation_source <- "User-confirmed"

        x <- rbind(x, newrow)

      } else {

        # If the confirmed compound is among candidates, elevate that row.
        h2 <- hit[x$compound_id[hit] == conf$compound_id[i]]

        if (!length(h2)) {
          # Preserve other candidates when adding a confirmed compound with a new ID.
          newrow <- x[hit[1],,drop=FALSE]
          newrow$compound_id <- conf$compound_id[i]
          newrow$metabolite_name <- conf$metabolite_name[i]
          newrow$annotation_score <- 1
          newrow$annotation_level <- "1_confirmed_standard"
          newrow$confirmed <- TRUE
          newrow$confirmation_source <- if ("source" %in% names(conf))
            as.character(conf$source[i]) else "User-confirmed"
          x <- rbind(x, newrow)
        } else {
          x$annotation_score[h2] <- 1
          x$annotation_level[h2] <- "1_confirmed_standard"
          x$confirmed[h2] <- TRUE
          x$confirmation_source[h2] <- if ("source" %in% names(conf))
            as.character(conf$source[i]) else "User-confirmed"
        }
      }
    }
  }

  # ---------------------------------------------------------------
  # Rank candidates per feature.
  # Ties are resolved by evidence hierarchy, score, ppm error, RT error.
  # ---------------------------------------------------------------
  level_rank <- c(
    "1_confirmed_standard" = 1,
    "2_msms_library" = 2,
    "3_mass_rt" = 3,
    "4_mass_only" = 4
  )

  x$level_rank <- unname(level_rank[x$annotation_level])

  x$abs_ppm_error <- abs(x$ppm_error)
  x$sort_rt_error <- ifelse(is.finite(x$rt_error), x$rt_error, Inf)

  ord <- order(
    x$mz_time,
    x$level_rank,
    -x$annotation_score,
    x$abs_ppm_error,
    x$sort_rt_error,
    na.last = TRUE
  )

  x <- x[ord,,drop=FALSE]

  x$candidate_rank <- ave(
    seq_len(nrow(x)),
    x$mz_time,
    FUN = seq_along
  )

  x$n_candidates <- ave(
    x$compound_id,
    x$mz_time,
    FUN = length
  )

  x$best_candidate <- x$candidate_rank == 1L

  # Flag features with multiple candidate rows, regardless of the top-ranked level.
  x$ambiguous <- x$n_candidates > 1

  rownames(x) <- NULL
  x
}



# -----------------------------------------------------------------------------
# ADDUCT AND EMPIRICAL-COMPOUND UTILITIES
# -----------------------------------------------------------------------------

pm_adduct_table <- function(mode=c("positive","negative")) {
  mode <- match.arg(mode)

  if (mode == "positive") {
    data.frame(
      adduct=c("[M+H]+","[M+Na]+","[M+NH4]+","[M+K]+","[M+2H]2+","[2M+H]+"),
      charge=c(1,1,1,1,2,1),
      multiplier=c(1,1,1,1,1,2),
      mass_shift=c(
        1.007276466621,
        22.989218,
        18.033823,
        38.963158,
        2*1.007276466621,
        1.007276466621
      ),
      stringsAsFactors=FALSE
    )
  } else {
    data.frame(
      adduct=c("[M-H]-","[M+Cl]-","[M+FA-H]-","[M+Ac-H]-","[M-2H]2-","[2M-H]-"),
      charge=c(-1,-1,-1,-1,-2,-1),
      multiplier=c(1,1,1,1,1,2),
      mass_shift=c(
        -1.007276466621,
        34.969402,
        44.998201,
        59.013851,
        -2*1.007276466621,
        -1.007276466621
      ),
      stringsAsFactors=FALSE
    )
  }
}

pm_neutral_mass_from_mz <- function(mz, adduct, mode=c("positive","negative")) {
  mode <- match.arg(mode)
  tab <- pm_adduct_table(mode)
  j <- match(adduct, tab$adduct)
  out <- rep(NA_real_, length(mz))

  ok <- !is.na(j) & is.finite(mz)
  if (any(ok)) {
    mult <- tab$multiplier[j[ok]]
    z <- abs(tab$charge[j[ok]])
    shift <- tab$mass_shift[j[ok]]
    out[ok] <- (mz[ok] * z - shift) / mult
  }
  out
}

pm_expand_adducts <- function(features,
                              mode=c("positive","negative"),
                              include_reported_adduct=TRUE) {
  mode <- match.arg(mode)
  f <- pm_read(features)

  if (!"mz_time" %in% names(f))
    stop("features must contain mz_time.")

  if (!"mz" %in% names(f))
    f$mz <- pm_parse_mz_time(f$mz_time)$mz

  tab <- pm_adduct_table(mode)

  out <- lapply(seq_len(nrow(f)), function(i) {
    z <- tab
    z$mz_time <- f$mz_time[i]
    # Repeat the scalar observed m/z once per adduct to preserve row alignment.
    z$observed_mz <- rep(as.numeric(f$mz[i]), nrow(z))

    z$neutral_mass <- pm_neutral_mass_from_mz(
      mz = z$observed_mz,
      adduct = z$adduct,
      mode = mode
    )

    if (length(z$neutral_mass) != nrow(z)) {
      stop(
        "Internal adduct-expansion error for ", f$mz_time[i],
        ": expected ", nrow(z), " neutral masses but obtained ",
        length(z$neutral_mass), "."
      )
    }

    z
  })

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL

  if (include_reported_adduct && "adduct" %in% names(f)) {
    rep_ad <- unique(f[,c("mz_time","adduct"),drop=FALSE])
    names(rep_ad)[2] <- "reported_adduct"
    ans <- merge(ans, rep_ad, by="mz_time", all.x=TRUE, sort=FALSE)
    ans$matches_reported_adduct <- ans$adduct == ans$reported_adduct
  }

  ans
}

pm_test_adduct_expansion <- function() {
  test <- data.frame(
    mz_time = c("72.01451_0.723"),
    mz = c(72.01451),
    adduct = c("[M-H]-"),
    stringsAsFactors = FALSE
  )

  neg <- pm_expand_adducts(
    test,
    mode = "negative",
    include_reported_adduct = TRUE
  )

  if (nrow(neg) != nrow(pm_adduct_table("negative")))
    stop("MetaEngine adduct self-test failed: wrong number of NEG hypotheses.")

  if (length(neg$neutral_mass) != nrow(neg))
    stop("MetaEngine adduct self-test failed: neutral-mass vector length mismatch.")

  if (any(!is.finite(neg$neutral_mass)))
    stop("MetaEngine adduct self-test failed: non-finite neutral mass.")

  invisible(TRUE)
}

pm_empirical_compounds_mztime <- function(X,
                                          feature_info=NULL,
                                          cor_thresh=0.80,
                                          rt_tol=0.15,
                                          method="spearman") {

  if (is.null(feature_info)) {
    feature_info <- pm_parse_mz_time(rownames(X))
  } else {
    feature_info <- pm_read(feature_info)
  }

  ids <- intersect(rownames(X), feature_info$mz_time)
  Xi <- X[ids,,drop=FALSE]
  fi <- feature_info[match(ids, feature_info$mz_time),,drop=FALSE]

  C <- suppressWarnings(
    cor(t(Xi), method=method, use="pairwise.complete.obs")
  )

  adj <- abs(C) >= cor_thresh
  rt <- as.numeric(fi$rt)
  rt_ok <- abs(outer(rt,rt,"-")) <= rt_tol
  adj <- adj & rt_ok
  diag(adj) <- TRUE

  if (requireNamespace("igraph", quietly=TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(
      adj, mode="undirected", diag=FALSE
    )
    cmp <- igraph::components(g)$membership
  } else {
    n <- nrow(adj)
    cmp <- integer(n)
    k <- 0
    for (i in seq_len(n)) {
      if (cmp[i] != 0) next
      k <- k + 1
      q <- i
      cmp[i] <- k
      while (length(q)) {
        u <- q[1]
        q <- q[-1]
        nb <- which(adj[u,] & cmp == 0)
        if (length(nb)) {
          cmp[nb] <- k
          q <- c(q, nb)
        }
      }
    }
  }

  data.frame(
    mz_time=ids,
    empirical_compound=paste0("EC",cmp),
    stringsAsFactors=FALSE
  )
}

pm_clean_metabolite_name <- function(x) {
  z <- as.character(x)
  z <- sub("^w/o MS2:", "", z, ignore.case=TRUE)
  z <- trimws(z)
  tolower(z)
}

pm_standardize_database <- function(db,
                                    database_name="database",
                                    saliva_database=FALSE) {

  x <- pm_read(db)

  # Flexible input aliases
  aliases <- list(
    compound_id=c("compound_id","accession","hmdb_id","HMDB_ID","LM_ID","id"),
    metabolite_name=c("metabolite_name","name","NAME","common_name"),
    exact_mass=c(
      "exact_mass","monoisotopic_molecular_weight",
      "monoisotopic_mass","EXACT_MASS","exact mass"
    ),
    formula=c("formula","chemical_formula","FORMULA"),
    inchikey=c("inchikey","inchi_key","INCHIKEY","INCHI_KEY"),
    smiles=c("smiles","SMILES"),
    kegg_id=c("kegg_id","KEGG_ID"),
    pubchem_id=c("pubchem_id","pubchem_cid","PUBCHEM_CID"),
    chebi_id=c("chebi_id","CHEBI_ID")
  )

  get_alias <- function(keys) {
    hit <- keys[keys %in% names(x)]
    if (length(hit)) return(hit[1])

    low <- tolower(names(x))
    for (k in keys) {
      j <- match(tolower(k), low)
      if (!is.na(j)) return(names(x)[j])
    }
    NULL
  }

  out <- data.frame(stringsAsFactors=FALSE)

  for (nm in names(aliases)) {
    src <- get_alias(aliases[[nm]])
    if (is.null(src)) {
      out[[nm]] <- if (nm == "exact_mass") rep(NA_real_, nrow(x)) else rep(NA_character_, nrow(x))
    } else {
      out[[nm]] <- x[[src]]
    }
  }

  out$compound_id <- as.character(out$compound_id)
  out$metabolite_name <- as.character(out$metabolite_name)
  out$exact_mass <- suppressWarnings(as.numeric(out$exact_mass))
  out$formula <- as.character(out$formula)
  out$inchikey <- as.character(out$inchikey)
  out$smiles <- as.character(out$smiles)

  if (all(is.na(out$compound_id)) && !all(is.na(out$metabolite_name)))
    out$compound_id <- out$metabolite_name

  out$database <- database_name
  out$saliva_database <- isTRUE(saliva_database)

  # Keep valid exact-mass records only
  out <- out[
    is.finite(out$exact_mass) &
    out$exact_mass > 0 &
    !is.na(out$metabolite_name) &
    nzchar(trimws(out$metabolite_name)),
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL
  out
}

pm_read_hmdb_xml <- function(xml_file,
                             database_name="HMDB",
                             saliva_database=TRUE) {

  pm_require("xml2")

  if (!file.exists(xml_file))
    stop("HMDB XML file not found: ", xml_file)

  doc <- xml2::read_xml(xml_file)

  # namespace-agnostic selection
  mets <- xml2::xml_find_all(doc, ".//*[local-name()='metabolite']")

  if (!length(mets))
    stop("No <metabolite> records found in HMDB XML.")

  val <- function(node, tag) {
    z <- xml2::xml_find_first(
      node,
      paste0("./*[local-name()='", tag, "']")
    )
    if (inherits(z, "xml_missing")) return(NA_character_)
    xml2::xml_text(z, trim=TRUE)
  }

  out <- data.frame(
    compound_id=vapply(mets, val, character(1), tag="accession"),
    metabolite_name=vapply(mets, val, character(1), tag="name"),
    exact_mass=suppressWarnings(as.numeric(
      vapply(mets, val, character(1), tag="monisotopic_molecular_weight")
    )),
    formula=vapply(mets, val, character(1), tag="chemical_formula"),
    inchikey=vapply(mets, val, character(1), tag="inchikey"),
    smiles=vapply(mets, val, character(1), tag="smiles"),
    kegg_id=NA_character_,
    pubchem_id=NA_character_,
    chebi_id=NA_character_,
    database=database_name,
    saliva_database=isTRUE(saliva_database),
    stringsAsFactors=FALSE
  )

  # HMDB has used both spellings in exports; repair exact mass if necessary.
  if (all(!is.finite(out$exact_mass))) {
    out$exact_mass <- suppressWarnings(as.numeric(
      vapply(mets, val, character(1), tag="monoisotopic_molecular_weight")
    ))
  }

  out <- out[
    is.finite(out$exact_mass) &
    out$exact_mass > 0,
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL
  out
}

pm_read_metabolite_database <- function(database_file,
                                        database_type=c("auto","csv","tsv","hmdb_xml"),
                                        database_name=NULL,
                                        saliva_database=FALSE) {

  database_type <- match.arg(database_type)

  if (is.data.frame(database_file)) {
    if (is.null(database_name)) database_name <- "local_database"
    return(
      pm_standardize_database(
        database_file,
        database_name=database_name,
        saliva_database=saliva_database
      )
    )
  }

  if (!file.exists(database_file))
    stop("Metabolite database file not found: ", database_file)

  ext <- tolower(tools::file_ext(database_file))

  if (database_type == "auto") {
    database_type <- if (ext == "xml") {
      "hmdb_xml"
    } else if (ext == "csv") {
      "csv"
    } else {
      "tsv"
    }
  }

  if (is.null(database_name)) {
    database_name <- if (database_type == "hmdb_xml") "HMDB" else "local_database"
  }

  if (database_type == "hmdb_xml") {
    return(
      pm_read_hmdb_xml(
        database_file,
        database_name=database_name,
        saliva_database=saliva_database
      )
    )
  }

  raw <- if (database_type == "csv") {
    read.csv(database_file, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    read.delim(database_file, check.names=FALSE, stringsAsFactors=FALSE)
  }

  pm_standardize_database(
    raw,
    database_name=database_name,
    saliva_database=saliva_database
  )
}

pm_match_database_by_neutral_mass <- function(features,
                                              database,
                                              mode=c("positive","negative"),
                                              ppm=5,
                                              reported_adduct=NULL) {

  mode <- match.arg(mode)
  f <- pm_read(features)

  if (!"mz_time" %in% names(f))
    stop("features require mz_time.")

  if (!"mz" %in% names(f))
    f$mz <- pm_parse_mz_time(f$mz_time)$mz

  if (!"rt" %in% names(f))
    f$rt <- pm_parse_mz_time(f$mz_time)$rt

  if (!is.null(reported_adduct)) {
    ra <- pm_read(reported_adduct)
    if (all(c("mz_time","adduct") %in% names(ra))) {
      f$adduct <- ra$adduct[match(f$mz_time, ra$mz_time)]
    }
  }

  adduct_hyp <- pm_expand_adducts(
    f,
    mode=mode,
    include_reported_adduct=("adduct" %in% names(f))
  )

  db <- database
  masses <- db$exact_mass

  out <- vector("list", nrow(adduct_hyp))

  for (i in seq_len(nrow(adduct_hyp))) {

    nmass <- adduct_hyp$neutral_mass[i]

    if (!is.finite(nmass) || nmass <= 0)
      next

    ppm_err <- 1e6 * (nmass - masses) / masses

    hit <- which(
      is.finite(ppm_err) &
      abs(ppm_err) <= ppm
    )

    if (!length(hit))
      next

    z <- db[hit,,drop=FALSE]

    z$mz_time <- adduct_hyp$mz_time[i]
    z$observed_mz <- adduct_hyp$observed_mz[i]
    z$observed_rt <- f$rt[
      match(adduct_hyp$mz_time[i], f$mz_time)
    ]
    z$adduct_hypothesis <- adduct_hyp$adduct[i]
    z$neutral_mass_hypothesis <- nmass
    z$database_ppm_error <- ppm_err[hit]

    if ("reported_adduct" %in% names(adduct_hyp)) {
      z$reported_adduct <- adduct_hyp$reported_adduct[i]
      z$matches_reported_adduct <- (
        !is.na(z$reported_adduct) &
        z$adduct_hypothesis == z$reported_adduct
      )
    } else {
      z$reported_adduct <- NA_character_
      z$matches_reported_adduct <- FALSE
    }

    out[[i]] <- z
  }

  ans <- do.call(rbind, out)

  if (is.null(ans)) {
    return(data.frame(
      mz_time=character(),
      compound_id=character(),
      metabolite_name=character(),
      exact_mass=numeric(),
      adduct_hypothesis=character(),
      database_ppm_error=numeric(),
      stringsAsFactors=FALSE
    ))
  }

  rownames(ans) <- NULL
  ans
}

pm_corroborate_database_candidates <- function(candidates,
                                               msdial_best=NULL,
                                               confirmed_annotations=NULL,
                                               ppm=5) {

  x <- candidates

  if (!nrow(x))
    return(x)

  x$msdial_preliminary_name <- NA_character_
  x$msdial_annotation_level <- NA_character_
  x$msdial_annotation_score <- NA_real_
  x$msdial_name_agreement <- FALSE
  x$msdial_inchikey_agreement <- FALSE
  x$msdial_formula_agreement <- FALSE
  x$msdial_msms_supported <- FALSE

  if (!is.null(msdial_best)) {

    m <- pm_read(msdial_best)
    idx <- match(x$mz_time, m$mz_time)

    if ("metabolite_name" %in% names(m))
      x$msdial_preliminary_name <- as.character(m$metabolite_name[idx])

    if ("annotation_level" %in% names(m))
      x$msdial_annotation_level <- as.character(m$annotation_level[idx])

    if ("annotation_score" %in% names(m))
      x$msdial_annotation_score <- as.numeric(m$annotation_score[idx])

    x$msdial_name_agreement <- (
      !is.na(x$msdial_preliminary_name) &
      pm_clean_metabolite_name(x$metabolite_name) ==
      pm_clean_metabolite_name(x$msdial_preliminary_name)
    )

    if ("inchikey" %in% names(m) && "inchikey" %in% names(x)) {
      mi <- as.character(m$inchikey[idx])
      x$msdial_inchikey_agreement <- (
        !is.na(x$inchikey) &
        !is.na(mi) &
        nzchar(x$inchikey) &
        x$inchikey == mi
      )
    }

    if ("formula" %in% names(m) && "formula" %in% names(x)) {
      mf <- as.character(m$formula[idx])
      x$msdial_formula_agreement <- (
        !is.na(x$formula) &
        !is.na(mf) &
        x$formula == mf
      )
    }

    if ("msms_evidence_score" %in% names(m)) {
      mm <- as.numeric(m$msms_evidence_score[idx])
      x$msdial_msms_supported <- is.finite(mm) & mm >= 0.60 &
        (
          x$msdial_name_agreement |
          x$msdial_inchikey_agreement
        )
    }
  }

  # Base database evidence from exact neutral-mass agreement.
  x$database_mass_score <- pmax(
    0,
    1 - abs(x$database_ppm_error) / ppm
  )

  x$reported_adduct_score <- ifelse(
    x$matches_reported_adduct %in% TRUE,
    1,
    0
  )

  x$saliva_context_score <- ifelse(
    "saliva_database" %in% names(x) & x$saliva_database %in% TRUE,
    1,
    0
  )

  # Candidate-ranking weights include a 0.05 contribution from saliva membership.
  x$database_annotation_score <- (
    0.60 * x$database_mass_score +
    0.10 * x$reported_adduct_score +
    0.05 * x$saliva_context_score +
    0.10 * ifelse(x$msdial_formula_agreement,1,0) +
    0.15 * ifelse(
      x$msdial_name_agreement | x$msdial_inchikey_agreement,
      1,
      0
    )
  )

  x$annotation_source <- "MetaEngine_database"
  x$annotation_level <- "4_database_mass_candidate"

  x$annotation_level[
    x$msdial_formula_agreement &
    !x$msdial_msms_supported
  ] <- "3_database_mass_plus_formula"

  x$annotation_level[
    (x$msdial_name_agreement | x$msdial_inchikey_agreement) &
    !x$msdial_msms_supported
  ] <- "3_database_MSDIAL_agreement"

  x$annotation_level[
    x$msdial_msms_supported
  ] <- "2_database_plus_msms"

  # Apply user-supplied standard confirmations.
  x$confirmed <- FALSE
  x$confirmation_source <- NA_character_

  conf <- pm_read_confirmed_annotations(confirmed_annotations)

  if (!is.null(conf)) {
    ci <- match(x$mz_time, conf$mz_time)
    hit <- which(!is.na(ci))

    if (length(hit)) {
      same <- pm_clean_metabolite_name(x$metabolite_name[hit]) ==
        pm_clean_metabolite_name(conf$metabolite_name[ci[hit]])

      h2 <- hit[same]

      if (length(h2)) {
        x$annotation_level[h2] <- "1_confirmed_standard"
        x$database_annotation_score[h2] <- 1
        x$confirmed[h2] <- TRUE
        x$confirmation_source[h2] <- if ("source" %in% names(conf))
          as.character(conf$source[ci[h2]]) else "User-confirmed"
      }
    }
  }

  # Rank independently derived candidates.
  level_rank <- c(
    "1_confirmed_standard"=1,
    "2_database_plus_msms"=2,
    "3_database_MSDIAL_agreement"=3,
    "3_database_mass_plus_formula"=4,
    "4_database_mass_candidate"=5
  )

  x$level_rank <- unname(level_rank[x$annotation_level])

  ord <- order(
    x$mz_time,
    x$level_rank,
    -x$database_annotation_score,
    abs(x$database_ppm_error),
    x$adduct_hypothesis,
    na.last=TRUE
  )

  x <- x[ord,,drop=FALSE]

  x$candidate_rank <- ave(
    seq_len(nrow(x)),
    x$mz_time,
    FUN=seq_along
  )

  x$n_database_candidates <- ave(
    x$compound_id,
    x$mz_time,
    FUN=length
  )

  x$best_candidate <- x$candidate_rank == 1L
  x$ambiguous <- x$n_database_candidates > 1L

  rownames(x) <- NULL
  x
}

pm_annotate_database <- function(features,
                                 database_file,
                                 database_type="auto",
                                 database_name=NULL,
                                 saliva_database=FALSE,
                                 ion_mode=c("positive","negative"),
                                 ppm=5,
                                 msdial_best=NULL,
                                 confirmed_annotations=NULL) {

  ion_mode <- match.arg(ion_mode)

  db <- pm_read_metabolite_database(
    database_file=database_file,
    database_type=database_type,
    database_name=database_name,
    saliva_database=saliva_database
  )

  message(
    "MetaEngine database annotation: loaded ",
    nrow(db), " compounds from ",
    if (nrow(db)) unique(db$database)[1] else "empty database", "."
  )

  candidates <- pm_match_database_by_neutral_mass(
    features=features,
    database=db,
    mode=ion_mode,
    ppm=ppm,
    reported_adduct=msdial_best
  )

  candidates <- pm_corroborate_database_candidates(
    candidates,
    msdial_best=msdial_best,
    confirmed_annotations=confirmed_annotations,
    ppm=ppm
  )

  if (!nrow(candidates)) {
    best <- data.frame(
      mz_time=as.character(features$mz_time),
      metabolite_name=NA_character_,
      annotation_source="MetaEngine_database",
      annotation_level="unannotated",
      database_annotation_score=NA_real_,
      n_database_candidates=0L,
      ambiguous=FALSE,
      stringsAsFactors=FALSE
    )
  } else {
    best <- candidates[
      candidates$best_candidate,
      ,
      drop=FALSE
    ]

    missing_ids <- setdiff(
      as.character(features$mz_time),
      best$mz_time
    )

    if (length(missing_ids)) {
      add <- data.frame(
        mz_time=missing_ids,
        metabolite_name=NA_character_,
        annotation_source="MetaEngine_database",
        annotation_level="unannotated",
        database_annotation_score=NA_real_,
        n_database_candidates=0L,
        ambiguous=FALSE,
        stringsAsFactors=FALSE
      )

      for (nm in setdiff(names(best), names(add)))
        add[[nm]] <- NA

      add <- add[,names(best),drop=FALSE]
      best <- rbind(best,add)
    }

    best <- best[
      match(as.character(features$mz_time), best$mz_time),
      ,
      drop=FALSE
    ]
  }

  summary <- data.frame(
    metric=c(
      "Total features",
      "Features with >=1 database candidate",
      "Features with multiple database candidates",
      "Features corroborated by MS-DIAL name/InChIKey",
      "Features supported by MS/MS"
    ),
    value=c(
      nrow(features),
      sum(best$n_database_candidates > 0, na.rm=TRUE),
      sum(best$n_database_candidates > 1, na.rm=TRUE),
      sum(
        (best$msdial_name_agreement %in% TRUE) |
        (best$msdial_inchikey_agreement %in% TRUE),
        na.rm=TRUE
      ),
      sum(best$annotation_level == "2_database_plus_msms", na.rm=TRUE)
    ),
    stringsAsFactors=FALSE
  )

  invisible(list(
    database=db,
    candidates=candidates,
    best=best,
    summary=summary
  ))
}

pm_write_database_annotation_outputs <- function(obj, output_dir) {

  write.csv(
    obj$candidates,
    file.path(output_dir,"database_annotation_all_candidates.csv"),
    row.names=FALSE
  )

  write.csv(
    obj$best,
    file.path(output_dir,"database_annotation_best_match.csv"),
    row.names=FALSE
  )

  write.csv(
    obj$summary,
    file.path(output_dir,"database_annotation_summary.csv"),
    row.names=FALSE
  )

  invisible(obj)
}

# -----------------------------------------------------------------------------
# MS-DIAL ANNOTATION IMPORT
# -----------------------------------------------------------------------------

pm_null_to_na <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    z <- trimws(x)
    z[tolower(z) %in% c("", "null", "na", "n/a", "nan", "none")] <- NA_character_
    return(z)
  }
  x
}

pm_num <- function(x) {
  x <- pm_null_to_na(x)
  suppressWarnings(as.numeric(x))
}

pm_flag <- function(x) {
  if (is.logical(x)) return(x)
  x <- pm_null_to_na(x)
  y <- tolower(as.character(x))
  out <- rep(NA, length(y))
  out[y %in% c("true","t","yes","y","1","matched")] <- TRUE
  out[y %in% c("false","f","no","n","0","unmatched")] <- FALSE
  out
}

pm_score01 <- function(x) {
  # Normalize common MS-DIAL score scales while preserving NA.
  # 0-1 stays 0-1; 0-100 becomes /100; 0-1000 becomes /1000.
  v <- pm_num(x)
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v)
  if (!any(ok)) return(out)

  vv <- v[ok]
  q99 <- suppressWarnings(as.numeric(stats::quantile(vv, 0.99, na.rm=TRUE, names=FALSE)))

  divisor <- if (is.finite(q99) && q99 <= 1.05) {
    1
  } else if (is.finite(q99) && q99 <= 105) {
    100
  } else {
    1000
  }

  out[ok] <- pmin(pmax(v[ok] / divisor, 0), 1)
  out
}

pm_has_spectrum <- function(x) {
  z <- pm_null_to_na(x)
  !is.na(z) & nzchar(trimws(as.character(z)))
}

pm_find_col <- function(dat, candidates, required=FALSE) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit)) return(hit[1])

  # Case-insensitive fallback
  low <- tolower(names(dat))
  for (cand in candidates) {
    j <- match(tolower(cand), low)
    if (!is.na(j)) return(names(dat)[j])
  }

  if (required)
    stop("Required MS-DIAL column not found. Tried: ",
         paste(candidates, collapse=", "))

  NULL
}

pm_read_msdial <- function(annotation_file,
                           sheet=NULL,
                           mz_time_digits=c(5,3)) {

  if (is.data.frame(annotation_file)) {
    raw <- annotation_file
    used_sheet <- NA_character_
  } else {
    if (!is.character(annotation_file) || length(annotation_file) != 1L)
      stop("annotation_file must be a data.frame or a file path.")

    if (!file.exists(annotation_file))
      stop("MS-DIAL annotation file not found: ", annotation_file)

    ext <- tolower(tools::file_ext(annotation_file))

    if (ext %in% c("xlsx","xls")) {
      pm_require("readxl")

      sheets <- readxl::excel_sheets(annotation_file)

      if (is.null(sheet)) {
        # Select the primary alignment or peak-height sheet.
        preferred <- grep(
          "height|alignment|aligned",
          sheets,
          ignore.case=TRUE,
          value=TRUE
        )

        preferred <- preferred[
          !grepl("stats|statistic", preferred, ignore.case=TRUE)
        ]

        used_sheet <- if (length(preferred)) preferred[1] else sheets[1]
      } else {
        used_sheet <- sheet
      }

      raw <- as.data.frame(
        readxl::read_excel(
          annotation_file,
          sheet=used_sheet,
          .name_repair="minimal"
        ),
        check.names=FALSE,
        stringsAsFactors=FALSE
      )

    } else if (ext == "csv") {
      raw <- read.csv(
        annotation_file,
        check.names=FALSE,
        stringsAsFactors=FALSE
      )
      used_sheet <- NA_character_

    } else {
      raw <- read.delim(
        annotation_file,
        check.names=FALSE,
        stringsAsFactors=FALSE
      )
      used_sheet <- NA_character_
    }
  }

  # Locate standard MS-DIAL columns.
  c_mztime <- pm_find_col(raw, c("mz_time","Mz_time","mz.time"))
  c_mz <- pm_find_col(raw, c("Average Mz","Average m/z","mz","m/z"))
  c_rt <- pm_find_col(raw, c("Average Rt(min)","Average RT(min)","time","rt","RT"))

  if (is.null(c_mztime)) {
    if (is.null(c_mz) || is.null(c_rt))
      stop("MS-DIAL annotation requires mz_time OR Average Mz + Average Rt(min).")

    mzv <- pm_num(raw[[c_mz]])
    rtv <- pm_num(raw[[c_rt]])

    raw$mz_time <- sprintf(
      paste0("%.", mz_time_digits[1], "f_%.", mz_time_digits[2], "f"),
      mzv, rtv
    )
    c_mztime <- "mz_time"
  }

  # Create standardized annotation table while retaining MS-DIAL fields.
  out <- data.frame(
    mz_time = as.character(raw[[c_mztime]]),
    stringsAsFactors=FALSE
  )

  out$mz <- if (!is.null(c_mz)) pm_num(raw[[c_mz]]) else pm_parse_mz_time(out$mz_time)$mz
  out$rt <- if (!is.null(c_rt)) pm_num(raw[[c_rt]]) else pm_parse_mz_time(out$mz_time)$rt

  pull <- function(cands, default=NA_character_) {
    nm <- pm_find_col(raw, cands)
    if (is.null(nm)) rep(default, nrow(raw)) else pm_null_to_na(raw[[nm]])
  }

  pull_num <- function(cands) {
    nm <- pm_find_col(raw, cands)
    if (is.null(nm)) rep(NA_real_, nrow(raw)) else pm_num(raw[[nm]])
  }

  pull_flag <- function(cands) {
    nm <- pm_find_col(raw, cands)
    if (is.null(nm)) rep(NA, nrow(raw)) else pm_flag(raw[[nm]])
  }

  out$metabolite_name <- as.character(pull(c("Metabolite name","Metabolite Name","Name")))
  out$adduct <- as.character(pull(c("Adduct type","Adduct Type","Adduct")))
  out$post_curation_result <- as.character(pull(c("Post curation result","Post curation")))
  out$msms_assigned <- pull_flag(c("MS/MS assigned","MSMS assigned"))
  out$reference_rt <- pull_num(c("Reference RT","Reference Rt"))
  out$reference_mz <- pull_num(c("Reference m/z","Reference Mz","Reference mz"))
  out$formula <- as.character(pull(c("Formula")))
  out$ontology <- as.character(pull(c("Ontology")))
  out$inchikey <- as.character(pull(c("INCHIKEY","InChIKey","Inchikey")))
  out$smiles <- as.character(pull(c("SMILES","Smiles")))
  out$annotation_tag <- as.character(pull(c("Annotation tag (VS1.0)","Annotation tag")))
  out$rt_matched <- pull_flag(c("RT matched","Rt matched"))
  out$mz_matched <- pull_flag(c("m/z matched","mz matched","M/Z matched"))
  out$msms_matched <- pull_flag(c("MS/MS matched","MSMS matched"))
  out$total_score_raw <- pull_num(c("Total score","Total Score"))
  out$rt_similarity_raw <- pull_num(c("RT similarity","RT Similarity"))
  out$dot_product_raw <- pull_num(c("Dot product","Dot Product"))
  out$reverse_dot_product_raw <- pull_num(c("Reverse dot product","Reverse Dot Product"))
  out$fragment_presence_raw <- pull_num(c("Fragment presence %","Fragment presence","Fragment Presence %"))
  out$sn_average <- pull_num(c("S/N average","S/N Average"))
  out$spectrum_reference_file <- as.character(pull(c("Spectrum reference file name","Spectrum reference file")))
  out$ms1_isotopic_spectrum <- as.character(pull(c("MS1 isotopic spectrum","MS1 Isotopic Spectrum")))
  out$msms_spectrum <- as.character(pull(c("MS/MS spectrum","MSMS spectrum")))

  # Normalize score-like evidence.
  out$total_score <- pm_score01(out$total_score_raw)
  out$rt_similarity <- pm_score01(out$rt_similarity_raw)
  out$dot_product <- pm_score01(out$dot_product_raw)
  out$reverse_dot_product <- pm_score01(out$reverse_dot_product_raw)
  out$fragment_presence <- pm_score01(out$fragment_presence_raw)

  out$has_ms1_isotope <- pm_has_spectrum(out$ms1_isotopic_spectrum)
  out$has_msms_spectrum <- pm_has_spectrum(out$msms_spectrum)

  # Treat generic Unknown names as unannotated.
  out$preliminary_id <- !is.na(out$metabolite_name) &
    !tolower(trimws(out$metabolite_name)) %in% c(
      "unknown","unk","unidentified","no match"
    )

  attr(out, "msdial_sheet") <- used_sheet
  out
}

pm_score_msdial_annotation <- function(msdial,
                                       confirmed_annotations=NULL,
                                       high_msms=0.80,
                                       moderate_msms=0.60) {

  x <- msdial

  # Encode positive evidence flags as 1 and missing flags as 0.
  mz_e <- ifelse(x$mz_matched %in% TRUE, 1, 0)

  # If Reference m/z exists, use ppm agreement as additional mass evidence.
  ppm_error <- rep(NA_real_, nrow(x))
  ok_mass <- is.finite(x$mz) & is.finite(x$reference_mz)
  ppm_error[ok_mass] <- pm_ppm(x$mz[ok_mass], x$reference_mz[ok_mass])

  ppm_score <- rep(0, nrow(x))
  ppm_score[ok_mass] <- pmax(0, 1 - abs(ppm_error[ok_mass]) / 10)
  mass_score <- pmax(mz_e, ppm_score, na.rm=TRUE)

  rt_flag <- ifelse(x$rt_matched %in% TRUE, 1, 0)
  rt_score <- pmax(
    rt_flag,
    ifelse(is.finite(x$rt_similarity), x$rt_similarity, 0),
    na.rm=TRUE
  )

  # Combine spectral scores, match flags, and assignment flags.
  msms_flag <- ifelse(x$msms_matched %in% TRUE, 1, 0)
  assigned_flag <- ifelse(x$msms_assigned %in% TRUE, 1, 0)

  spectral_components <- cbind(
    ifelse(is.finite(x$dot_product), x$dot_product, 0),
    ifelse(is.finite(x$reverse_dot_product), x$reverse_dot_product, 0),
    ifelse(is.finite(x$fragment_presence), x$fragment_presence, 0)
  )

  weighted_spectral <- (
    0.45 * spectral_components[,1] +
    0.35 * spectral_components[,2] +
    0.20 * spectral_components[,3]
  )

  msms_score <- pmax(
    weighted_spectral,
    ifelse(is.finite(x$total_score), x$total_score, 0),
    0.90 * msms_flag,
    0.30 * assigned_flag,
    na.rm=TRUE
  )

  isotope_score <- ifelse(x$has_ms1_isotope, 1, 0)

  # The weighted evidence score does not establish an MSI identification level.
  x$annotation_score <- (
    0.30 * mass_score +
    0.15 * rt_score +
    0.45 * msms_score +
    0.10 * isotope_score
  )

  x$ppm_error <- ppm_error
  x$mass_evidence_score <- mass_score
  x$rt_evidence_score <- rt_score
  x$msms_evidence_score <- msms_score
  x$isotope_evidence_score <- isotope_score

  # Evidence-derived confidence tiers.
  x$annotation_level <- "unannotated"

  # Preliminary exact-mass name only.
  x$annotation_level[x$preliminary_id & mass_score > 0] <- "4_mass_candidate"

  # Accurate mass + RT.
  x$annotation_level[
    x$preliminary_id &
    mass_score > 0 &
    rt_score >= 0.50
  ] <- "3_mass_rt"

  # Assign MS/MS evidence tiers using score thresholds and available spectrum flags.
  x$annotation_level[
    x$preliminary_id &
    mass_score > 0 &
    msms_score >= moderate_msms
  ] <- "2B_msms_supported"

  x$annotation_level[
    x$preliminary_id &
    mass_score > 0 &
    msms_score >= high_msms &
    (x$msms_matched %in% TRUE | x$has_msms_spectrum)
  ] <- "2A_high_confidence_msms"

  x$confirmed <- FALSE
  x$confirmation_source <- NA_character_

  # Authentic/internal-standard confirmation overrides computational evidence.
  conf <- pm_read_confirmed_annotations(confirmed_annotations)

  if (!is.null(conf)) {
    idx <- match(x$mz_time, conf$mz_time)
    hit <- which(!is.na(idx))

    if (length(hit)) {
      x$metabolite_name[hit] <- conf$metabolite_name[idx[hit]]
      if ("compound_id" %in% names(conf))
        x$compound_id <- if ("compound_id" %in% names(x)) x$compound_id else NA_character_

      if ("compound_id" %in% names(conf)) {
        if (!"compound_id" %in% names(x))
          x$compound_id <- NA_character_
        x$compound_id[hit] <- conf$compound_id[idx[hit]]
      }

      x$annotation_level[hit] <- "1_confirmed_standard"
      x$annotation_score[hit] <- 1
      x$confirmed[hit] <- TRUE
      x$confirmation_source[hit] <- if ("source" %in% names(conf))
        as.character(conf$source[idx[hit]]) else "User-confirmed"
    }
  }

  # Preserve one MS-DIAL row per mz_time. If duplicates exist, choose strongest
  # evidence but retain duplicate count.
  x$n_msdial_rows <- ave(x$mz_time, x$mz_time, FUN=length)

  level_rank <- c(
    "1_confirmed_standard"=1,
    "2A_high_confidence_msms"=2,
    "2B_msms_supported"=3,
    "3_mass_rt"=4,
    "4_mass_candidate"=5,
    "unannotated"=6
  )

  x$level_rank <- unname(level_rank[x$annotation_level])

  ord <- order(
    x$mz_time,
    x$level_rank,
    -x$annotation_score,
    na.last=TRUE
  )
  x <- x[ord,,drop=FALSE]
  x$best_candidate <- !duplicated(x$mz_time)

  rownames(x) <- NULL
  x
}

pm_annotate_msdial <- function(annotation_file,
                               confirmed_annotations=NULL,
                               sheet=NULL,
                               high_msms=0.80,
                               moderate_msms=0.60) {

  imported <- pm_read_msdial(
    annotation_file=annotation_file,
    sheet=sheet
  )

  scored <- pm_score_msdial_annotation(
    imported,
    confirmed_annotations=confirmed_annotations,
    high_msms=high_msms,
    moderate_msms=moderate_msms
  )

  best <- scored[scored$best_candidate,,drop=FALSE]

  level_order <- c(
    "1_confirmed_standard",
    "2A_high_confidence_msms",
    "2B_msms_supported",
    "3_mass_rt",
    "4_mass_candidate",
    "unannotated"
  )

  tab <- table(factor(best$annotation_level, levels=level_order))

  summary <- data.frame(
    annotation_level=level_order,
    n_features=as.integer(tab),
    percent=round(100*as.integer(tab)/max(nrow(best),1),2),
    stringsAsFactors=FALSE
  )

  summary <- rbind(
    data.frame(
      annotation_level="Total_features",
      n_features=nrow(best),
      percent=100,
      stringsAsFactors=FALSE
    ),
    data.frame(
      annotation_level="Any_preliminary_ID",
      n_features=sum(best$preliminary_id,na.rm=TRUE),
      percent=round(100*mean(best$preliminary_id,na.rm=TRUE),2),
      stringsAsFactors=FALSE
    ),
    data.frame(
      annotation_level="MSMS_assigned",
      n_features=sum(best$msms_assigned %in% TRUE,na.rm=TRUE),
      percent=round(100*mean(best$msms_assigned %in% TRUE,na.rm=TRUE),2),
      stringsAsFactors=FALSE
    ),
    data.frame(
      annotation_level="MSMS_matched",
      n_features=sum(best$msms_matched %in% TRUE,na.rm=TRUE),
      percent=round(100*mean(best$msms_matched %in% TRUE,na.rm=TRUE),2),
      stringsAsFactors=FALSE
    ),
    summary
  )

  invisible(list(
    source="MSDIAL",
    imported=imported,
    candidates=scored,
    best=best,
    summary=summary,
    sheet=attr(imported,"msdial_sheet")
  ))
}

pm_annotate <- function(features,
                        compound_library,
                        confirmed_annotations = NULL,
                        msms_evidence = NULL,
                        ppm = 5,
                        rt_tol = 0.30,
                        msms_score_threshold = 0.70,
                        mz_time_col = "mz_time",
                        compound_id_col = "compound_id",
                        compound_name_col = "metabolite_name",
                        compound_mz_col = "mz",
                        compound_rt_col = "rt") {

  candidates <- pm_match_mass(
    features = features,
    compounds = compound_library,
    ppm = ppm,
    rt_tol = rt_tol,
    mz_time_col = mz_time_col,
    compound_id_col = compound_id_col,
    compound_name_col = compound_name_col,
    compound_mz = compound_mz_col,
    compound_rt = compound_rt_col
  )

  # Include confirmed annotations when the mass library returns no candidates.
  conf <- pm_read_confirmed_annotations(confirmed_annotations)

  if (!nrow(candidates) && !is.null(conf)) {
    shell <- data.frame(
      mz_time = conf$mz_time,
      compound_id = conf$compound_id,
      metabolite_name = conf$metabolite_name,
      observed_mz = pm_parse_mz_time(conf$mz_time)$mz,
      observed_rt = pm_parse_mz_time(conf$mz_time)$rt,
      library_mz = NA_real_,
      library_rt = NA_real_,
      ppm_error = NA_real_,
      rt_error = NA_real_,
      rt_matched = FALSE,
      stringsAsFactors = FALSE
    )
    candidates <- shell
  }

  scored <- pm_score_annotations(
    candidates,
    confirmed_annotations = confirmed_annotations,
    msms_evidence = msms_evidence,
    ppm = ppm,
    rt_tol = rt_tol,
    msms_score_threshold = msms_score_threshold
  )

  all_features <- pm_annotation_features(features, mz_time_col=mz_time_col)

  if (is.null(scored) || !nrow(scored)) {
    best <- all_features
    best$compound_id <- NA_character_
    best$metabolite_name <- NA_character_
    best$annotation_level <- "unannotated"
    best$annotation_score <- NA_real_
    best$n_candidates <- 0L
    best$ambiguous <- FALSE
    best$best_candidate <- TRUE
  } else {

    best <- scored[scored$best_candidate,,drop=FALSE]

    # Add features with no annotation candidate.
    missing_ids <- setdiff(all_features$mz_time, best$mz_time)

    if (length(missing_ids)) {
      u <- all_features[match(missing_ids, all_features$mz_time),,drop=FALSE]

      add <- data.frame(
        mz_time = u$mz_time,
        compound_id = NA_character_,
        metabolite_name = NA_character_,
        annotation_level = "unannotated",
        annotation_score = NA_real_,
        n_candidates = 0L,
        ambiguous = FALSE,
        best_candidate = TRUE,
        stringsAsFactors = FALSE
      )

      # Align columns with best.
      for (nm in setdiff(names(best), names(add)))
        add[[nm]] <- NA

      add <- add[,names(best),drop=FALSE]
      best <- rbind(best, add)
    }

    best <- best[match(all_features$mz_time, best$mz_time),,drop=FALSE]
  }

  # Summary by feature, not candidate row.
  level_order <- c(
    "1_confirmed_standard",
    "2_msms_library",
    "3_mass_rt",
    "4_mass_only",
    "unannotated"
  )

  counts <- table(
    factor(best$annotation_level, levels=level_order)
  )

  summary <- data.frame(
    annotation_level = level_order,
    n_features = as.integer(counts),
    percent = round(
      100 * as.integer(counts) / nrow(best),
      2
    ),
    stringsAsFactors=FALSE
  )

  summary <- rbind(
    data.frame(
      annotation_level="Total_features",
      n_features=nrow(best),
      percent=100,
      stringsAsFactors=FALSE
    ),
    data.frame(
      annotation_level="Any_annotation",
      n_features=sum(best$annotation_level != "unannotated", na.rm=TRUE),
      percent=round(
        100 * mean(best$annotation_level != "unannotated", na.rm=TRUE),
        2
      ),
      stringsAsFactors=FALSE
    ),
    data.frame(
      annotation_level="Ambiguous_best_annotations",
      n_features=sum(best$ambiguous, na.rm=TRUE),
      percent=round(
        100 * mean(best$ambiguous, na.rm=TRUE),
        2
      ),
      stringsAsFactors=FALSE
    ),
    summary
  )

  invisible(
    list(
      candidates = scored,
      best = best,
      summary = summary
    )
  )
}

pm_merge_annotation_results <- function(stats, annotation_best) {

  st <- pm_read(stats)
  an <- pm_read(annotation_best)

  if (!"mz_time" %in% names(st) && "feature_id" %in% names(st))
    names(st)[names(st)=="feature_id"] <- "mz_time"

  if (!"mz_time" %in% names(st))
    stop("Differential results require mz_time.")

  if (!"mz_time" %in% names(an))
    stop("Annotation results require mz_time.")

  idx <- match(st$mz_time, an$mz_time)

  preferred <- c(
    "metabolite_name","compound_id","adduct","formula","ontology","inchikey",
    "smiles","annotation_tag","annotation_level","annotation_score",
    "confirmed","confirmation_source","preliminary_id",
    "reference_mz","reference_rt","ppm_error",
    "mz_matched","rt_matched","msms_assigned","msms_matched",
    "total_score_raw","rt_similarity_raw","dot_product_raw",
    "reverse_dot_product_raw","fragment_presence_raw",
    "mass_evidence_score","rt_evidence_score","msms_evidence_score",
    "isotope_evidence_score","has_ms1_isotope","has_msms_spectrum",
    "spectrum_reference_file","ms1_isotopic_spectrum","msms_spectrum",
    "n_candidates","ambiguous","n_msdial_rows"
  )

  keep_cols <- intersect(preferred, names(an))

  for (nm in keep_cols)
    st[[nm]] <- an[[nm]][idx]

  st
}

pm_write_annotation_outputs <- function(annotation_object, output_dir) {

  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

  cand <- annotation_object$candidates
  best <- annotation_object$best
  smry <- annotation_object$summary

  write.csv(
    cand,
    file.path(output_dir, "annotation_all_candidates.csv"),
    row.names=FALSE
  )

  write.csv(
    best,
    file.path(output_dir, "annotation_best_match.csv"),
    row.names=FALSE
  )

  write.csv(
    smry,
    file.path(output_dir, "annotation_summary.csv"),
    row.names=FALSE
  )

  if (!is.null(annotation_object$source) &&
      identical(annotation_object$source, "MSDIAL")) {

    evidence_cols <- intersect(
      c(
        "mz_time","metabolite_name","adduct","formula","ontology","inchikey",
        "annotation_tag","reference_mz","reference_rt",
        "mz_matched","rt_matched","msms_assigned","msms_matched",
        "total_score_raw","rt_similarity_raw","dot_product_raw",
        "reverse_dot_product_raw","fragment_presence_raw",
        "ms1_isotopic_spectrum","msms_spectrum",
        "mass_evidence_score","rt_evidence_score","msms_evidence_score",
        "isotope_evidence_score","annotation_score","annotation_level",
        "confirmed","confirmation_source"
      ),
      names(best)
    )

    write.csv(
      best[,evidence_cols,drop=FALSE],
      file.path(output_dir, "annotation_evidence_MSDIAL.csv"),
      row.names=FALSE
    )
  }

  lines <- c(
    "============================================================",
    "                 MetaEngine Annotation Summary",
    "============================================================"
  )

  if (!is.null(annotation_object$source))
    lines <- c(lines, paste0("Annotation source: ", annotation_object$source))

  if (!is.null(annotation_object$sheet) &&
      !is.na(annotation_object$sheet))
    lines <- c(lines, paste0("Imported sheet: ", annotation_object$sheet))

  lines <- c(lines, "")

  for (i in seq_len(nrow(smry))) {
    lines <- c(
      lines,
      sprintf(
        "%-32s %7d (%6.2f%%)",
        smry$annotation_level[i],
        smry$n_features[i],
        smry$percent[i]
      )
    )
  }

  writeLines(
    lines,
    file.path(output_dir, "annotation_summary.txt")
  )

  invisible(annotation_object)
}




