

test_that("parse.seq works correctly",{

  #basic test
  x <- parse.seq("GPPGPPGPPGPKGPPGPPGPPGPP")
  expect_match(x$seq[1], "GPPGPPGPPGPK")

  input <- "GPPGLPGPPGPPKPPPPPPFPPPRPPPPPPPPPPPPVPPRKRQQQQQQQQQQQQQQQQQIQR"

  #cut before
  psb <- parse.seq(input,cutbefore=T,verbose=T)
  expect_match(psb$seq[1], "GPPGLPGPPGPP")

  psb <- parse.seq(input,cuts = "L|F|V|I|A",cutbefore = T)
  psb <- unique(psb$seq)
  expect_match(psb[1], "LPGPPGPPKPPPPPP")
  expect_match(psb[2], "FPPPRPPPPPPPPPPPP")


})
