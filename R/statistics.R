# MetaEngine statistics module.
# model.matrix() encodes continuous and categorical covariates.
# block_col uses limma::duplicateCorrelation for repeated or paired observations.

pm_differential <- function(X, metadata,
                            formula,
                            contrast=NULL,
                            block_col=NULL,
                            robust=TRUE,
                            trend=FALSE) {
  pm_require("limma")
  if (!inherits(formula, "formula")) formula <- as.formula(formula)
  
  design <- model.matrix(formula, data=metadata)
  
  original_design_names <- colnames(design)
  colnames(design) <- make.names(colnames(design), unique=TRUE)
  
  if (!identical(original_design_names, colnames(design))) {
    changed <- which(original_design_names != colnames(design))
    if (length(changed)) {
      message("MetaEngine: sanitized design coefficient names:")
      for (i in changed)
        message("  ", original_design_names[i], " -> ", colnames(design)[i])
    }
  }
  block <- NULL
  corfit <- NULL
  if (!is.null(block_col)) {
    if (!block_col %in% names(metadata)) stop("block_col not found.")
    block <- metadata[[block_col]]
    corfit <- limma::duplicateCorrelation(X, design=design, block=block)
    fit <- limma::lmFit(X, design, block=block, correlation=corfit$consensus)
  } else {
    fit <- limma::lmFit(X, design)
  }
  
  if (!is.null(contrast)) {
    C <- if (is.character(contrast))
      limma::makeContrasts(contrasts=contrast, levels=design) else contrast
    fit <- limma::contrasts.fit(fit, C)
  }
  fit <- limma::eBayes(fit, robust=robust, trend=trend)
  tab <- limma::topTable(fit, number=Inf, sort.by="none", adjust.method="BH")
  tab$mz_time <- rownames(tab)
  tab <- tab[,c("mz_time",setdiff(names(tab),"mz_time")),drop=FALSE]
  
  attr(tab,"design") <- design
  attr(tab,"duplicateCorrelation") <- if (!is.null(corfit)) corfit$consensus else NA_real_
  tab
}

# Mixed-model analysis using variancePartition::dream.

pm_differential_dream <- function(X, metadata, formula, contrast=NULL) {
  pm_require("variancePartition")
  pm_require("limma")
  if (!inherits(formula,"formula")) formula <- as.formula(formula)
  vobj <- limma::voom(X, model.matrix(~1, metadata), plot=FALSE)
  fit <- variancePartition::dream(vobj, formula, metadata)
  if (!is.null(contrast)) {
    L <- variancePartition::makeContrastsDream(formula, metadata, contrasts=contrast)
    fit <- variancePartition::dream(vobj, formula, metadata, L=L)
  }
  fit <- limma::eBayes(fit)
  tab <- limma::topTable(fit, number=Inf, sort.by="none")
  tab$feature_id <- rownames(tab)
  tab[,c("feature_id",setdiff(names(tab),"feature_id")),drop=FALSE]
}


