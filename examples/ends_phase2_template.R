# Set the local input file paths before running this example.
library(MetaEngine)
intensity_file <- "path/to/intensity_mz_time.tsv"
metadata_file <- "path/to/metadata.csv"
msdial_file <- "path/to/MSDIAL_annotation.xlsx"
hmdb_file <- "path/to/hmdb_metabolites.xml"
meta <- read.csv(metadata_file, check.names=FALSE)
meta$Exposure <- relevel(factor(meta$Exposure), ref="Nonvaper")
for (nm in intersect(c("Sex","Race","BMI_ind","Sugar_ind","EDU_ind","OHDNTVST_ind"), names(meta))) meta[[nm]] <- factor(meta[[nm]])
args <- list(
  intensity=intensity_file, metadata=meta, contrast="ExposureVaper",
  feature_id_col="mz_time", sample_id_col="sample_id", group_col="Exposure",
  annotation_source="MSDIAL", annotation_file=msdial_file,
  metabolite_database=hmdb_file, metabolite_database_type="hmdb_xml",
  metabolite_database_name="HMDB", ion_mode="negative", database_ppm=5,
  expand_adducts=TRUE, empirical_compounds=TRUE,
  max_missing=0.5, max_group_missing=0.8, impute=TRUE, log2=TRUE, quantile=TRUE
)
unadjusted <- do.call(pm_run, c(args, list(formula=~Exposure, output_dir="outputs/unadjusted")))
adjusted <- do.call(pm_run, c(args, list(formula=~Exposure+Age+Sex+BMI_ind+Race+Sugar_ind+EDU_ind+OHDNTVST_ind, output_dir="outputs/adjusted")))
# Each model runs annotation for its retained features.
head(adjusted$stats[, intersect(c("mz_time","P.Value","adj.P.Val","metabolite_name","annotation_level"), names(adjusted$stats))])
