[README.md](https://github.com/user-attachments/files/35164241/README.md)
# MetaEngine

MetaEngine is an R package for post-acquisition metabolomics. It uses `mz_time` as the feature identifier and provides accurate-mass database candidate generation alongside MS-DIAL preliminary identification evidence.

## Modules

| File | Responsibility | Status |
| --- | --- | --- |
| `R/annotation.R` | HMDB XML/CSV import, ion/adduct matching, MS-DIAL evidence, candidate ranking | Active development |
| `R/annotation_unified.R` | Unified annotation fields in main results | Development |
| `R/knowledge_base.R` | Optional saliva/oral/exposure context | Development; requires curated evidence |
| `R/statistics.R` | Covariate-adjusted differential models | Experimental, independently developed |
| `R/pathways.R` | ORA, rank scoring, overlap/stability diagnostics | Experimental, independently developed |
| `R/preprocessing.R` | Post-acquisition feature filtering, imputation, transforms | Experimental |
| `R/workflow.R` | Orchestration via `pm_run()` | Experimental |

No raw peak picking or chromatographic alignment is provided.

## Installation

Install R >= 4.2.0. Clone this repository, then in R:

```r
install.packages("remotes")
remotes::install_local("/path/to/MetaEngine", dependencies = TRUE)
library(MetaEngine)
```

`limma` and `preprocessCore` are Bioconductor packages. Install separately if needed:

```r
install.packages("BiocManager")
BiocManager::install(c("limma", "preprocessCore"), ask = FALSE, update = FALSE)
```

`pm_initialize()` checks dependencies. `pm_initialize(install = TRUE)` installs missing dependencies.

## Database annotation example

```r
library(MetaEngine)
features <- data.frame(mz_time = c("181.07066_5.000", "150.00000_2.000"))
db <- data.frame(compound_id = "TEST001", metabolite_name = "Synthetic glucose-like compound",
                 exact_mass = 180.063383, formula = "C6H12O6")
res <- pm_annotate_database(
  features = pm_parse_mz_time(features$mz_time),
  database_file = db, ion_mode = "positive", ppm = 10
)
head(res$candidates)
```

The example uses a synthetic reference record to demonstrate mass matching.

HMDB XML files are supplied separately and can be passed to `pm_annotate_database()` or `pm_run()`. Source terms and provenance requirements are described in [Data sources](docs/DATA_SOURCES.md).

## Unified result fields

`pm_run()` retains `msdial_*` and `metaengine_*` provenance fields and adds:

`mz_time`, `metabolite_name`, `annotation_source`, `annotation_level`, `annotation_score`.

Mass matches represent candidate annotations. Annotation scores are heuristic rankings, not identification probabilities. Saliva occurrence and other contextual associations do not establish analytical identity or increase identification confidence.

The existing `pm_run()` argument `feature_id_col = "mz_time"` is retained for compatibility.

## Analysis examples

- `examples/annotation_demo.R`: database candidate search using synthetic data.
- `examples/ends_phase2_template.R`: cohort workflow with configurable local input paths.

## Development and tests

```r
install.packages("testthat")
testthat::test_local("/path/to/MetaEngine")
```

GitHub Actions is configured to run R package checks on pushes and pull requests. A full `R CMD check` and end-to-end validation remain pending. Statistics and pathway analysis are developed and validated as separate components. See [Known limitations](docs/KNOWN_LIMITATIONS.md).

## Citation and licensing

The MIT license applies to MetaEngine code. HMDB, MS-DIAL, and other databases and software retain their own terms. See [Data sources](docs/DATA_SOURCES.md).
