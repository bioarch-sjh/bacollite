?

test_that("GPM structure is valid",{

  #Load gpm from text file
  f<-NA

  #Load saved gpm data
  s<-NA

  expect_equal(nrow(s),nrow(f))
  expect_equal(ncol(s),ncol(f))


})
