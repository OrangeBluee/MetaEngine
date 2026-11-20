# MetaEngine pathways module.

pm_pathway_feature_table <- function(feature_map,pathway_map) {
  fm <- pm_prepare_feature_mapping(feature_map)
  pm <- pm_read(pathway_map)
  if (!all(c("compound_id","pathway_id") %in% names(pm)))
    stop("pathway_map needs compound_id and pathway_id.")
  z <- merge(fm,pm,by="compound_id",all=FALSE)
  if (!"pathway_name" %in% names(z)) z$pathway_name <- z$pathway_id
  z
}

# -----------------------------------------------------------------------------
# PATHWAY OVERREPRESENTATION ANALYSIS
# Feature-level resampling from a fixed background of mapped features.
# -----------------------------------------------------------------------------

pm_pathway_ora <- function(stats, feature_map, pathway_map,
                           p_cut=0.05, p_col="P.Value",
                           permutations=2000, seed=1,
                           min_hits=3) {
  st <- pm_read(stats)
  if (!all(c("feature_id",p_col) %in% names(st))) stop("stats columns missing.")
  pft <- pm_pathway_feature_table(feature_map,pathway_map)
  
  universe <- intersect(st$feature_id,unique(pft$feature_id))
  sig <- intersect(st$feature_id[is.finite(st[[p_col]]) & st[[p_col]]<p_cut],universe)
  paths <- unique(pft[,c("pathway_id","pathway_name")])
  
  score_one <- function(selected) {
    vapply(paths$pathway_id,function(pid){
      members <- unique(pft$feature_id[pft$pathway_id==pid])
      a <- length(intersect(selected,members))
      b <- length(selected)-a
      c <- length(setdiff(intersect(universe,members),selected))
      d <- length(universe)-a-b-c
      if (a<1) return(0)
      -log10(fisher.test(matrix(c(a,b,c,d),2),alternative="greater")$p.value)
    },numeric(1))
  }
  
  obs <- score_one(sig)
  set.seed(seed)
  null <- matrix(0,nrow=permutations,ncol=length(obs))
  if (length(sig)>0) for (b in seq_len(permutations))
    null[b,] <- score_one(sample(universe,length(sig),replace=FALSE))
  
  emp <- vapply(seq_along(obs),function(j)
    (1+sum(null[,j]>=obs[j]))/(permutations+1),numeric(1))
  
  hits <- vapply(paths$pathway_id,function(pid)
    length(intersect(sig,unique(pft$feature_id[pft$pathway_id==pid]))),integer(1))
  
  out <- data.frame(paths,hits=hits,score=obs,p_empirical=emp)
  out$FDR <- pm_bh(out$p_empirical)
  out$reportable <- out$hits>=min_hits
  out[order(out$p_empirical,-out$hits),]
}

# -----------------------------------------------------------------------------
# WEIGHTED PATHWAY SCORING AND OVERLAP DIAGNOSTICS
# Score pathways using mapped feature evidence and normalized candidate weights.
# Related functions report pathway overlap, shared compounds, and cutoff sensitivity.
# -----------------------------------------------------------------------------

pm_rank_score <- function(stats,p_col="P.Value",effect_col=NULL) {
  st <- pm_read(stats)
  p <- pmax(as.numeric(st[[p_col]]), .Machine$double.xmin)
  s <- -log10(p)
  if (!is.null(effect_col) && effect_col %in% names(st))
    s <- s*sign(as.numeric(st[[effect_col]]))
  data.frame(feature_id=st$feature_id,rank_score=s)
}

pm_weighted_pathway_score <- function(stats,feature_map,pathway_map,
                                      p_col="P.Value",effect_col=NULL,
                                      permutations=5000,seed=1,
                                      min_features=3) {
  rs <- pm_rank_score(stats,p_col,effect_col)
  pft <- pm_pathway_feature_table(feature_map,pathway_map)
  z <- merge(pft,rs,by="feature_id")
  paths <- unique(z[,c("pathway_id","pathway_name")])
  
  # Calculate each pathway's weighted mean absolute feature score.
  calc <- function(vals) {
    vapply(paths$pathway_id,function(pid){
      q <- vals[z$pathway_id==pid]
      w <- z$map_weight[z$pathway_id==pid]
      if (length(unique(z$feature_id[z$pathway_id==pid]))<min_features) return(NA_real_)
      sum(w*abs(q),na.rm=TRUE)/sum(w[is.finite(q)],na.rm=TRUE)
    },numeric(1))
  }
  
  obs <- calc(z$rank_score)
  set.seed(seed)
  null <- matrix(NA_real_,nrow=permutations,ncol=length(obs))
  base <- rs$rank_score
  names(base) <- rs$feature_id
  for (b in seq_len(permutations)) {
    perm <- sample(base,length(base),replace=FALSE)
    names(perm) <- names(base)
    null[b,] <- calc(perm[z$feature_id])
  }
  p <- vapply(seq_along(obs),function(j) {
    if (!is.finite(obs[j])) return(NA_real_)
    v <- null[,j]; v<-v[is.finite(v)]
    (1+sum(v>=obs[j]))/(1+length(v))
  },numeric(1))
  
  nf <- vapply(paths$pathway_id,function(pid)
    length(unique(z$feature_id[z$pathway_id==pid])),integer(1))
  nc <- vapply(paths$pathway_id,function(pid)
    length(unique(z$compound_id[z$pathway_id==pid])),integer(1))
  
  out <- data.frame(paths,n_features=nf,n_candidate_compounds=nc,
                    stable_score=obs,p_stable=p,FDR_stable=pm_bh(p))
  out[order(out$p_stable,-out$stable_score),]
}

pm_pathway_overlap <- function(feature_map,pathway_map,metric=c("jaccard","overlap")) {
  metric <- match.arg(metric)
  z <- pm_pathway_feature_table(feature_map,pathway_map)
  ids <- unique(z$pathway_id)
  M <- matrix(0,length(ids),length(ids),dimnames=list(ids,ids))
  sets <- lapply(ids,function(p) unique(z$compound_id[z$pathway_id==p]))
  names(sets)<-ids
  for(i in seq_along(ids)) for(j in i:length(ids)) {
    a<-sets[[i]]; b<-sets[[j]]
    val <- if(metric=="jaccard") length(intersect(a,b))/length(union(a,b)) else
      length(intersect(a,b))/min(length(a),length(b))
    if(!is.finite(val)) val<-0
    M[i,j]<-M[j,i]<-val
  }
  M
}

pm_pathway_families <- function(results,feature_map,pathway_map,
                                overlap_threshold=0.5) {
  R <- results
  O <- pm_pathway_overlap(feature_map,pathway_map,"jaccard")
  ids <- intersect(R$pathway_id,rownames(O))
  A <- O[ids,ids,drop=FALSE] >= overlap_threshold
  diag(A)<-TRUE
  if(requireNamespace("igraph",quietly=TRUE)){
    g<-igraph::graph_from_adjacency_matrix(A,mode="undirected",diag=FALSE)
    fam<-igraph::components(g)$membership
  } else {
    n<-nrow(A); fam<-integer(n);k<-0
    for(i in seq_len(n)) if(fam[i]==0){k<-k+1;q<-i;fam[i]<-k
    while(length(q)){u<-q[1];q<-q[-1];nb<-which(A[u,]&fam==0)
    if(length(nb)){fam[nb]<-k;q<-c(q,nb)}}}
  }
  map <- data.frame(pathway_id=ids,pathway_family=paste0("PF",fam))
  out <- merge(R,map,by="pathway_id",all.x=TRUE)
  
  # Select the family representative by p_stable, or p_empirical when available.
  pnm <- if("p_stable"%in%names(out)) "p_stable" else
    if("p_empirical"%in%names(out)) "p_empirical" else NULL
  out$family_representative <- FALSE
  if(!is.null(pnm)){
    for(f in unique(na.omit(out$pathway_family))){
      ix<-which(out$pathway_family==f)
      out$family_representative[ix[which.min(out[[pnm]][ix])]]<-TRUE
    }
  }
  out
}

pm_pathway_evidence <- function(results,feature_map,pathway_map) {
  z <- pm_pathway_feature_table(feature_map,pathway_map)
  selected <- results$pathway_id
  z <- z[z$pathway_id %in% selected,,drop=FALSE]
  
  # Count reported pathways supported by each mapped compound.
  mult <- aggregate(pathway_id~compound_id,z,function(x)length(unique(x)))
  names(mult)[2]<-"n_reported_pathways"
  z <- merge(z,mult,by="compound_id",all.x=TRUE)
  
  ev <- lapply(selected,function(pid){
    q<-z[z$pathway_id==pid,,drop=FALSE]
    compounds<-unique(q$compound_id)
    nmulti<-mult$n_reported_pathways[match(compounds,mult$compound_id)]
    data.frame(
      pathway_id=pid,
      n_compounds=length(compounds),
      n_unique_compounds=sum(nmulti==1,na.rm=TRUE),
      n_shared_compounds=sum(nmulti>1,na.rm=TRUE),
      uniqueness_fraction=mean(nmulti==1,na.rm=TRUE),
      shared_evidence_fraction=mean(nmulti>1,na.rm=TRUE)
    )
  })
  do.call(rbind,ev)
}

# Repeat ORA across feature p-value cutoffs using the same background universe.

pm_threshold_stability <- function(stats,feature_map,pathway_map,
                                   thresholds=c(0.001,0.005,0.01,0.05,0.1),
                                   p_col="P.Value",
                                   permutations=1000,seed=1,min_hits=3) {
  ans <- lapply(seq_along(thresholds),function(i){
    x<-pm_pathway_ora(stats,feature_map,pathway_map,p_cut=thresholds[i],
                      p_col=p_col,permutations=permutations,
                      seed=seed+i,min_hits=min_hits)
    x$feature_p_threshold<-thresholds[i]
    x
  })
  do.call(rbind,ans)
}

# Compare feature sets selected by nominal and adjusted p-value cutoffs.
# Report their sizes, differences, and whether the adjusted set is a subset.

pm_compare_feature_thresholds <- function(stats,p_cut=0.05,fdr_cut=0.05,
                                          p_col="P.Value",fdr_col="adj.P.Val") {
  st<-pm_read(stats)
  A<-st$feature_id[is.finite(st[[p_col]]) & st[[p_col]]<p_cut]
  B<-st$feature_id[is.finite(st[[fdr_col]]) & st[[fdr_col]]<fdr_cut]
  data.frame(
    nominal_n=length(A), fdr_n=length(B),
    fdr_not_in_nominal=length(setdiff(B,A)),
    nominal_not_in_fdr=length(setdiff(A,B)),
    fdr_subset_nominal=all(B%in%A)
  )
}


# -----------------------------------------------------------------------------
# FEATURE-TO-COMPOUND MAPPING
# feature_map columns: feature_id, compound_id, optional annotation_score.
# pathway_map columns: compound_id, pathway_id, optional pathway_name.
# Multiple candidate compounds per feature are retained.
# -----------------------------------------------------------------------------


pm_prepare_feature_mapping <- function(feature_map) {
  fm <- pm_read(feature_map)
  if (!all(c("feature_id","compound_id") %in% names(fm)))
    stop("feature_map needs feature_id and compound_id.")
  
  # Normalize candidate weights within each feature; use equal weights when
  # annotation_score is absent.
  if ("annotation_score" %in% names(fm)) {
    fm$raw_weight <- pmax(as.numeric(fm$annotation_score),1e-6)
  } else {
    fm$raw_weight <- 1
  }
  den <- ave(fm$raw_weight,fm$feature_id,FUN=sum)
  fm$map_weight <- fm$raw_weight/den
  fm
}

