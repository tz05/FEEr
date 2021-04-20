## ----setup--------------------------------------------------------------------
library(FEE)

## ---- include = FALSE---------------------------------------------------------
set.seed(100)
d_trait <- cbind(tr1 = rnorm(10,0,1), tr2 = rnorm(10,4,3))
d_comm <- rbind(comm1 = sample(0:5, 10, replace = TRUE, prob = 6:1),
                comm2 = sample(0:5, 10, replace = TRUE, prob = 6:1),
                comm3 = sample(0:5, 10, replace = TRUE, prob = 6:1),
                comm4 = sample(0:5, 10, replace = TRUE))
d_trait
d_comm

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

