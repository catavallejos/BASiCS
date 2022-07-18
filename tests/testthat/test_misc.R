test_that("BASiCS_EffectiveSize works", {
  data(ChainSC)
  ess <- BASiCS_EffectiveSize(ChainSC, Param = "mu")
  expect_type(ess, "double")
  expect_equal(ess[[1]], 129.73861)
  expect_warning(
    BASiCS_effectiveSize(ChainSC, Param = "mu"),
    "'BASiCS_effectiveSize' is deprecated."
  )
})

test_that("utils work", {
  data(ChainSC)
  expect_equal(BASiCS:::.DistanceName("ResDisp"), "distance")
  expect_equal(BASiCS:::.DistanceName("Mean"), "fold change")
  expect_equal(BASiCS:::.LogDistanceName("ResDisp"), "distance")
  expect_equal(BASiCS:::.LogDistanceName("Mean"), "log2(fold change)")
  expect_equal(BASiCS:::.cap("test"), "Test")
  expect_equal(BASiCS:::.NSamples(ChainSC), 75)
  expect_equal(BASiCS:::.MeasureName("Mean"), "mean expression")
  expect_equal(BASiCS:::.MeasureName("Disp"), "over dispersion")
  expect_equal(BASiCS:::.MeasureName("ResDisp"), "residual over dispersion")
})
