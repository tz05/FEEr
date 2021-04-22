## ---- message = FALSE------------------------------------------------------------------------------------------
library(FEE)

## --------------------------------------------------------------------------------------------------------------
set.seed(321)
d_trait <- cbind(tr1 = rnorm(10, 0, 1), tr2 = rnorm(10, 3, 2))
d_comm <- rbind(comm1 = sample(0:4, 10, replace = TRUE, prob = c(8, 4:1)),
                comm2 = sample(0:4, 10, replace = TRUE, prob = c(8, 4:1)),
                comm3 = sample(0:4, 10, replace = TRUE, prob = c(8, 4:1)),
                comm4 = sample(0:4, 10, replace = TRUE),
                comm5 = sample(0:4, 10, replace = TRUE),
                comm6 = sample(0:4, 10, replace = TRUE))
d_trait
d_comm

## ---- fig.height = 2.4, fig.width = 3.6------------------------------------------------------------------------
par(mfrow = c(2, 3), mar = rep(0.3, 4))
for(i in 1:6) {
  plot(d_trait, axes = F, type = "n")
  polygon(rbind(c(min(d_trait[,1]), min(d_trait[,2])), c(max(d_trait[,1]), min(d_trait[,2])),
                c(max(d_trait[,1]), max(d_trait[,2])), c(min(d_trait[,1]), max(d_trait[,2]))), lty=2)
  points(d_trait[which(d_comm[i,] > 0),], cex = 0.4*d_comm[i, which(d_comm[i,] > 0)], col = 2, pch = 19)
  legend("topright", legend = paste0("comm", i, " "), bty = "n")
}

## ---- results = "hide", message = FALSE------------------------------------------------------------------------
fee1 <- computeFEE(d_trait, d_comm)

## --------------------------------------------------------------------------------------------------------------
fee1

## ---- fig.width = 2, fig.height = 2----------------------------------------------------------------------------
nsp <- computeNSP(d_comm)
nsp
par(mar = c(2.8, 2.8, 0.2, 0.2), mgp = c(1.8,0.5,0))
plot(nsp, fee1, ylim = 0:1, yaxs = "i")

## ---- results = "hide", message = FALSE------------------------------------------------------------------------
fee2 <- computeFEE(d_trait, d_comm, abundWeighted = FALSE)
fee3 <- computeFEE(d_trait, d_comm, dis_metric = "manhattan")
fee4 <- computeFEE(d_trait, d_comm, poolBased = TRUE)
fee5 <- computeFEE(d_trait, d_comm, abundWeighted = FALSE, dis_metric = "manhattan")
fee6 <- computeFEE(d_trait, d_comm, abundWeighted = FALSE, poolBased = TRUE)
fee7 <- computeFEE(d_trait, d_comm, dis_metric = "manhattan", poolBased = TRUE)
fee8 <- computeFEE(d_trait, d_comm, abundWeighted = FALSE, dis_metric = "manhattan", poolBased = TRUE)

## ---- fig.width = 6, fig.height = 6----------------------------------------------------------------------------
fees <- cbind(fee1, fee2, fee3, fee4, fee5, fee6, fee7, fee8) # FEE under different parameters
panel_ln <- function(x,y) {   # show the 1:1 line in the scatter-plots
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(-0.02, 1.02, -0.02, 1.02))
  abline(0, 1, col = 2)
  points(x, y)
}
pairs(fees, gap = 0.3, panel = panel_ln)

## ---- results = "hide", message = FALSE------------------------------------------------------------------------
# calculate both FEE0 and FEE for comm2 and comm3
fee1b <- computeFEE(d_trait, d_comm[2:3,], doFEE0 = TRUE)

## --------------------------------------------------------------------------------------------------------------
fee1b

## ---- results = "hide", message = FALSE------------------------------------------------------------------------
if(file.exists("ecdf_exp")) file.remove("ecdf_exp")
str_fee4b <- "computeFEE(d_trait, d_comm, poolBased = TRUE, user_ecdf = \"ecdf_exp\")"
str_fee6b <- "computeFEE(d_trait, d_comm, abundWeighted = FALSE, poolBased = TRUE, user_ecdf = \"ecdf_exp\")"
t_fee4b <- system.time(eval(parse(text = str_fee4b)))   # calculate eCDFs from scratch and save them
t_fee6b <- system.time(eval(parse(text = str_fee6b)))   # use saved eCDFs

## --------------------------------------------------------------------------------------------------------------
rbind(summary(t_fee4b), summary(t_fee6b))

## --------------------------------------------------------------------------------------------------------------
# MST branch lengths of comm1, with and without species abundance considered
br1a <- communityMST(d_comm[1,], d_trait)
br1b <- communityMST(d_comm[1,], d_trait, abundWeighted = FALSE)
rbind(br1a,br1b)

## ---- results = "hide", message = FALSE------------------------------------------------------------------------
# build eCDF for the case of 4 species in 2D trait space
ecdf_sp4a <- FEE0_eCDF(4, 2)                            # no prior knowledge of trait values
ecdf_sp4b <- FEE0_eCDF(4, 2, pool_traits = d_trait)     # species pool is specified

## ---- fig.width = 4, fig.height = 2----------------------------------------------------------------------------
# comm2 and comm3 have 4 species, so use them to illustrate the calculation
fee0 <- fee1b$FEE0    # FEE0 of comm2 and comm3, see "fee1b" above
cbind(FEE0 = fee0, FEEa = ecdf_sp4a(fee0), FEEb = ecdf_sp4b(fee0))
par(mfrow = c(1,2), mar = c(2.8, 2.8, 0.2, 0.2), mgp = c(1.8,0.5,0))
plot(ecdf_sp4a, xlab = "FEE0", ylab = "FEEa", main = "", pch = ".")
points(fee0, ecdf_sp4a(fee0), pch = 19)
plot(ecdf_sp4b, xlab = "FEE0", ylab = "FEEb", main = "", pch = ".")
points(fee0, ecdf_sp4b(fee0), pch = 19)
fee1b

