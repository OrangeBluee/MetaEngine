test_that("mz_time parsing preserves IDs", {
  x <- pm_parse_mz_time(c("181.07066_5.000", "150.00000_2.000"))
  expect_identical(x$mz_time[1], "181.07066_5.000")
  expect_equal(x$mz[1], 181.07066)
  expect_equal(x$rt[2], 2)
})

test_that("adduct expansion yields one neutral mass per hypothesis", {
  f <- data.frame(mz_time="181.07066_5.000", mz=181.07066)
  z <- pm_expand_adducts(f, mode="positive")
  expect_equal(nrow(z), nrow(pm_adduct_table("positive")))
  expect_equal(length(z$neutral_mass), nrow(z))
})

test_that("database annotation generates candidates independently of MS-DIAL", {
  f <- pm_parse_mz_time("181.07066_5.000")
  db <- data.frame(compound_id="SYNTH001", metabolite_name="Synthetic mass standard", exact_mass=180.063383)
  z <- pm_annotate_database(f, db, ion_mode="positive", ppm=10)
  expect_gt(nrow(z$candidates), 0)
  expect_true("mz_time" %in% names(z$best))
})

test_that("unified output retains provenance and consistent main columns", {
  x <- data.frame(mz_time=c("a","b","c"), metaengine_metabolite_name=c("A",NA,NA),
                  msdial_metabolite_name=c("B","C",NA),
                  metaengine_annotation_level=c("4_database_mass_candidate",NA,NA),
                  msdial_evidence_level=c("2B_msms_supported","3_mass_rt",NA))
  y <- pm_unify_annotation_results(x)
  expect_identical(y$metabolite_name, c("A","C",NA_character_))
  expect_identical(y$annotation_source[1:2], c("MetaEngine_database","MSDIAL"))
  expect_true("msdial_metabolite_name" %in% names(y))
})
