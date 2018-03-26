

test_that("parse.seq works correctly",{

  #basic test
  x <- parse.seq("GPPGPPGPPGPKGPPGPPGPPGPP")
  expect_match(x$seq[1], "GPPGPPGPPGPK")

  input <- "GPPGPPGPPGPPKPPPPPPPPPPRPPPPPPPPPPPPPPPRKRQQQQQQQQQQQQQQQQQQQR"

  #cut before
  psb <- parse.seq(input,cutbefore=T,verbose=T)
  expect_match(psb$seq[1], "GPPGPPGPPGPP")



})
