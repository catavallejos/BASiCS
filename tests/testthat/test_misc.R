context("miscellaneous")

test_that("BASiCS_EffectiveSize works", {
  data(ChainSC)
  ess <- BASiCS_EffectiveSize(ChainSC, Param = "mu")
  expect_is(ess, "numeric")
  expect_equal(ess[[1]], 129.73861)
  expect_warning(
    BASiCS_effectiveSize(ChainSC, Param = "mu"),
    "'BASiCS_effectiveSize' is deprecated."
  )
})

test_that("utils work", {
  data(ChainSC)
  expect_equal(DistanceName("ResDisp"), "distance")
  expect_equal(DistanceName("Mean"), "fold change")
  expect_equal(LogDistanceName("ResDisp"), "distance")
  expect_equal(LogDistanceName("Mean"), "log2(fold change)")
  expect_equal(cap("test"), "Test")
  expect_equal(NSamples(ChainSC), 75)
  expect_equal(MeasureName("Mean"), "mean expression")
  expect_equal(MeasureName("Disp"), "over dispersion")
  expect_equal(MeasureName("ResDisp"), "residual over dispersion")
})
