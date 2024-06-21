library(psych)
library(rbenchmark)
library(microbenchmark)
X <- matrix(rnorm(10000), 100, 100)

mbm <- microbenchmark("lm" = { y <- tr(X%*%X) },
                      "pseudoinverse" = {
                        y <- tr(X%*%X%*%X)
                      },
                      times=1000)
mbm

