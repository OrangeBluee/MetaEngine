# MetaEngine preprocessing module.

pm_load_data <- function(intensity, metadata,
                         feature_id_col=1,
                         sample_id_col="sample_id",
                         annotation_cols=NULL) {
  dat <- pm_read(intensity)
  meta <- pm_read(metadata)
  
  fid_name <- if (is.numeric(feature_id_col)) names(dat)[feature_id_col] else feature_id_col
  if (!fid_name %in% names(dat)) stop("Feature ID column not found.")
  if (!sample_id_col %in% names(meta)) stop("sample_id column not found.")
  
  ids <- as.character(dat[[fid_name]])
  samp <- as.character(meta[[sample_id_col]])
  sample_cols <- intersect(samp, names(dat))
  if (length(sample_cols) < 2) stop("No/too few matching sample columns.")
  
  X <- as.matrix(dat[, sample_cols, drop=FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- ids
  meta <- meta[match(colnames(X), meta[[sample_id_col]]),,drop=FALSE]
  rownames(meta) <- meta[[sample_id_col]]
  
  if (is.null(annotation_cols)) {
    annotation_cols <- setdiff(names(dat), c(fid_name, sample_cols))
  }
  ann <- data.frame(feature_id=ids, stringsAsFactors=FALSE)
  if (length(annotation_cols)) ann <- cbind(ann, dat[, intersect(annotation_cols,names(dat)),drop=FALSE])
  
  list(X=X, metadata=meta, annotation=ann)
}

# -----------------------------------------------------------------------------
# PREPROCESSING
# -----------------------------------------------------------------------------

pm_filter_missing <- function(X, max_missing=0.5, group=NULL, max_group_missing=NULL) {
  keep <- rowMeans(is.na(X) | !is.finite(X) | X==0) <= max_missing
  if (!is.null(group) && !is.null(max_group_missing)) {
    g <- as.factor(group)
    kg <- sapply(levels(g), function(z) {
      rowMeans(is.na(X[,g==z,drop=FALSE]) | !is.finite(X[,g==z,drop=FALSE]) |
                 X[,g==z,drop=FALSE]==0) <= max_group_missing
    })
    if (is.vector(kg)) kg <- matrix(kg,ncol=1)
    keep <- keep & apply(kg,1,all)
  }
  X[keep,,drop=FALSE]
}

pm_impute_halfmin <- function(X) {
  Y <- X
  for (i in seq_len(nrow(Y))) {
    v <- Y[i,]
    pos <- v[is.finite(v) & v>0]
    replacement <- if (length(pos)) min(pos)/2 else NA_real_
    bad <- !is.finite(v) | is.na(v) | v==0
    Y[i,bad] <- replacement
  }
  Y
}

pm_transform <- function(X, log2=TRUE, quantile=FALSE, center=FALSE, scale=FALSE) {
  Y <- X
  if (log2) {
    if (any(Y <= 0, na.rm=TRUE)) stop("Non-positive values remain before log2.")
    Y <- log2(Y)
  }
  if (quantile) {
    pm_require("preprocessCore")
    rn <- rownames(Y); cn <- colnames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    rownames(Y)<-rn; colnames(Y)<-cn
  }
  if (center || scale) {
    Y <- t(scale(t(Y), center=center, scale=scale))
  }
  Y
}


