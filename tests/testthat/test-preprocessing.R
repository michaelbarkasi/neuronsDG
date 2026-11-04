# Regression values captured by running these functions against the real
# kilosort4 demo data (downloaded/cached per R/demo_data.R) on the current
# package build.

test_that("load_mat() reads a real includeVector.mat file (regression)", {
  mat_path <- file.path(
    kilo4_demo_path(),
    "2025-04-05_11-50-43", "kilosort4", "includeVector.mat"
  )
  out <- neuronsDG:::load_mat(mat_path)

  expect_equal(dim(out), c(102, 1))
  expect_true(all(out %in% c(0, 1)))
  expect_equal(sum(out), 76)
  expect_equal(as.vector(out)[1:10], c(1, 1, 1, 0, 0, 0, 1, 0, 1, 1))
})

test_that("import.kilo4() returns spikes/stim data frames with expected shape (regression)", {
  folder <- file.path(kilo4_demo_path(), "2025-04-05_11-50-43")
  res <- import.kilo4(
    folder_path = folder,
    stim_data_file = file.path(folder, "StimulusStamps.csv"),
    trial_time_start = 500,
    trial_time_end = 500 + 1520,
    verbose = FALSE
  )

  expect_named(res, c("spikes", "stim"))
  expect_equal(
    colnames(res$spikes),
    c("spike_clusters", "spike_times", "cluster_group", "stim_responsive",
      "trial_number", "trial_time")
  )
  expect_equal(nrow(res$spikes), 570256)
  expect_equal(dim(res$stim), c(500, 5))
  expect_equal(sum(!is.na(res$spikes$trial_number)), 207320)
})

test_that("preprocess.kilo4() reproduces the exact quality-control counts described in the kilosort tutorial", {
  # vignettes/tutorial_tau_est_kilosort.Rmd's narrative text states: "for
  # Mecp2 HET mice, 26 clusters from the right hemisphere of females passed
  # ... For Shank3 KO mice, 9 clusters passed from the left hemisphere of
  # males and 1 from the left hemisphere of females ... For C57 (wildtype)
  # mice, 1 cluster passed from the left hemisphere of females ... As the
  # last neuron, number 37, is the only C57 sample ..."
  processed <- kilo4_demo_preprocessed()

  expect_named(processed, c("spikes", "timeXtrial", "cluster.key"))
  expect_equal(nrow(processed$cluster.key), 37)

  summary <- summarize.cluster.key(
    key = processed$cluster.key,
    covariate_list = c("genotype", "hemi", "sex")
  )
  non_na <- summary[!is.na(summary$n_cells), ]
  non_na <- non_na[order(non_na$genotype, non_na$hemi, non_na$sex), ]

  expect_equal(nrow(non_na), 4)
  expect_equal(sum(non_na$n_cells), 37)

  mecp2 <- non_na[non_na$genotype == "Mecp2 HET", ]
  expect_equal(as.character(mecp2$hemi), "RH")
  expect_equal(as.character(mecp2$sex), "F")
  expect_equal(mecp2$n_cells, 26)

  shank3 <- non_na[non_na$genotype == "Shank3 KO", ]
  expect_equal(nrow(shank3), 2)
  expect_true(all(as.character(shank3$hemi) == "LH"))
  expect_setequal(shank3$n_cells, c(9, 1))

  c57 <- non_na[non_na$genotype == "C57", ]
  expect_equal(as.character(c57$hemi), "LH")
  expect_equal(as.character(c57$sex), "F")
  expect_equal(c57$n_cells, 1)
})

test_that("summarize.cluster.key() computes rounded mean spikes/trials per combination", {
  processed <- kilo4_demo_preprocessed()
  summary <- summarize.cluster.key(
    key = processed$cluster.key,
    covariate_list = c("genotype", "hemi", "sex")
  )
  mecp2 <- summary[summary$genotype == "Mecp2 HET" & summary$hemi == "RH" & summary$sex == "F", ]

  expect_equal(mecp2$n_cells, 26)
  expect_equal(mecp2$mean_spikes, round(mean(
    processed$cluster.key$num.of.spikes[processed$cluster.key$genotype == "Mecp2 HET"]
  ), 1))
})
