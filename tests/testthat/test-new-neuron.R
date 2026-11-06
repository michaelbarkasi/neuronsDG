test_that("new.neuron() uses its documented defaults", {
  nrn <- new.neuron()
  id <- fetch.id.data(nrn)

  expect_s3_class(nrn, "neuron")
  expect_equal(id$id_num, 0)
  expect_equal(id$recording_name, "not_provided")
  expect_equal(id$type, "generic")
  expect_equal(id$sex, "not_provided")
  expect_equal(id$hemi, "not_provided")
  expect_equal(id$region, "not_provided")
  expect_equal(id$age, "not_provided")
  expect_false(id$sim)
  expect_equal(id$unit_time, "ms")
  expect_equal(id$unit_sample_rate, "Hz")
  expect_equal(id$unit_data, "mV")
  expect_equal(id$t_per_bin, 10.0)
  expect_equal(id$sample_rate, 1e4)
})

test_that("new.neuron() stores custom arguments", {
  nrn <- new.neuron(
    id_num = 7,
    recording_name = "rec_A",
    type = "LIF",
    genotype = "KO",
    sex = "M",
    hemi = "right",
    region = "CA1",
    age = "P14",
    sim = TRUE,
    unit_time = "s",
    unit_sample_rate = "kHz",
    unit_data = "spike",
    t_per_bin = 5.0,
    sample_rate = 2e3
  )
  id <- fetch.id.data(nrn)

  expect_equal(id$id_num, 7)
  expect_equal(id$recording_name, "rec_A")
  expect_equal(id$type, "LIF")
  expect_equal(id$genotype, "KO")
  expect_equal(id$sex, "M")
  expect_equal(id$hemi, "right")
  expect_equal(id$region, "CA1")
  expect_equal(id$age, "P14")
  expect_true(id$sim)
  expect_equal(id$unit_time, "s")
  expect_equal(id$unit_sample_rate, "kHz")
  expect_equal(id$unit_data, "spike")
  expect_equal(id$t_per_bin, 5.0)
  expect_equal(id$sample_rate, 2e3)
})

test_that("each new.neuron() call is an independent object", {
  n1 <- new.neuron(id_num = 1)
  n2 <- new.neuron(id_num = 2)
  expect_equal(fetch.id.data(n1)$id_num, 1)
  expect_equal(fetch.id.data(n2)$id_num, 2)
})
