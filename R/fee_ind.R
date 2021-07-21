globalVariables(c('coor_grids','i'))

#' Calculate Functional Extension and Evenness (FEE) index
#'
#' Computes the FEE index for the communities specified in \code{abund}, based
#' on the trait data specified in \code{pool_traits}.
#' 
#' This function is the core function of the package. The \code{pool_traits}
#' and \code{abund} are key arguments. The \code{dis_metric},
#' \code{abundWeighted}, and \code{poolBased} are essential arguments, which
#' together define the FEE calculation scheme. The other three arguments are
#' auxiliary arguments, which might affect calculation efficiency or alter
#' output format, but do not affect FEE calculation results.
#' 
#' Preparing eCDFs usually takes the most time of FEE calculation. To increase
#' calculation speed, some pre-calculated eCDFs are enclosed in the package.
#' For the cases with 1 to 6 traits, 2 to 90 species, and under the assumption
#' of no prior knowledge about trait distribution (i.e., \code{poolBased} is
#' \code{FALSE}, thus species pool does not affect eCDFs), these pre-calculated
#' eCDFs are used in the FEE calculation. For other scenarios, especially
#' species-pool-based calculation (\code{poolBased = TRUE}), eCDFs are not
#' prepared and need to be calculated from scratch. In this case, specifying
#' \code{user_ecdf} allows the newly calculated eCDFs to be saved and reused
#' as needed. As a result, the saved eCDFs can considerably reduce calculation
#' time for communities that are under the same settings of \code{pool_traits},
#' \code{dis_metric}, and \code{poolBased}.
#' 
#' @param pool_traits a matrix (or data frame) of numeric functional traits for
#'   species pool. Columns are for different traits and rows are for different
#'   species in the species pool.
#' @param abund matrix (or data frame) of species abundances in the communities
#'   of interest. Columns are for different species in the species pool (in the
#'   same species order as in \code{pool_traits}. Thus its number of columns
#'   needs to be equal to the number of rows in \code{pool_traits}). Each row
#'   is for a community. For species not presented in a community, their
#'   abundance values should be assigned 0 or \code{NA}.
#' @param dis_metric string specifying the scheme of quantifying species
#'   dissimilarity in trait space. The currently available options are
#'   \code{"euclidean"} (the default) and \code{"manhattan"}. Euclidean
#'   distances are root sum-of-squares of differences, and manhattan distances
#'   are the sum of absolute differences.
#' @param abundWeighted logical; indicating whether the abundance values in
#'   \code{abund} are considered as weights in FEE calculation. The default is
#'   \code{TRUE}.
#' @param poolBased logical; indicating whether the calculation is based on
#'   species pool. The default is \code{FALSE}, indicating no prior knowledge
#'   about trait distribution (i.e., adopting uniform distribution for trait
#'   values to build null communities). If \code{TRUE}, \code{pool_traits}
#'   is used to build null communities.
#' @param doFEE0 logical; indicating whether the output includes FEE0 (raw
#'   functional extension and evenness metric). The default is \code{FALSE}.
#' @param user_ecdf string specifying the path and name of the file that saves
#'   user's results of empirical cumulative distribution function (eCDF) of
#'   FEE0. Under the same setting (\code{dis_metric} and \code{pool_traits}),
#'   calculated eCDFs could be reused to save time. If the argument is
#'   \code{NULL} (default), no user's eCDF will be adopted or saved, and the
#'   species-pool-based eCDFs will be calculated from scratch. If the file path
#'   is provided, the eCDFs contained in this file will be loaded and used as
#'   needed, and when the eCDF of a specific combination of species number and
#'   trait number is not available in the file, it will be calculated from
#'   scratch and then saved to the file.
#' @param doParallel logical; indicating whether the calculation of FEE is
#'   conducted in parallel. The default value is \code{FALSE}, i.e.,
#'   sequentially calculating the FEE index for the communities specified in
#'   \code{abund}.
#' @return The FEE index values of the communities. If \code{doFEE0} is
#'   \code{TRUE}, the output includes two lists: the raw FEE0 metrics and the
#'   FEE indices. Specifically, FEE values for monocultures are \code{NA}.
#' @examples
#' d_trait <- cbind(tr1 = rnorm(10, 0, 1), tr2 = rnorm(10, 3, 2))
#' d_comm <- rbind(comm1 = c(0:9),
#'                 comm2 = rep(1, 10),
#'                 comm3 = runif(10, 0, 5),
#'                 comm4 = rchisq(10, 2))
#' computeFEE(d_trait, d_comm, abundWeighted = TRUE)
#' @import doSNOW
#' @import foreach
#' @export
computeFEE <- function(pool_traits, abund, dis_metric = c("euclidean","manhattan"), abundWeighted = TRUE,
                       poolBased = FALSE, doFEE0 = FALSE, user_ecdf = NULL, doParallel = FALSE) {
  dis_metric <- match.arg(dis_metric)
  # Assertion statements on input
  if(is.vector(abund)) {
    stopifnot(nrow(pool_traits) == length(abund))
    abund <- matrix(abund,nrow=1)
  } else stopifnot(nrow(pool_traits) == ncol(abund))
  stopifnot(dis_metric %in% c("euclidean","manhattan"))
  
  num_traits <- ncol(pool_traits)
  abund[is.na(abund)] <- 0
  nsp_comm <- computeNSP(abund)
  n_core_slurm <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
  numCores <- parallel::detectCores()
  progress <- function(n) utils::setTxtProgressBar(pb,n)
  
  std_pool_traits <- apply(pool_traits,2,function(x) (x-min(x))/(max(x)-min(x)))  # standardize traits to 0-1
  message(sprintf("\nCalculating MSTs for the %d communities...",length(nsp_comm)))
  pb <- utils::txtProgressBar(max = nrow(abund), initial = NA, style = 3)
  if(doParallel & nrow(abund)>1) {
    cl <- snow::makeCluster(min(numCores,n_core_slurm,nrow(abund),na.rm=T))
    registerDoSNOW(cl)
    message(sprintf("With %d cores...",min(numCores,n_core_slurm,nrow(abund),na.rm=T)))
    progress(0)
    #registerDoParallel()
    opts <- list(progress = progress)
    comm_MST <- foreach(i=1:nrow(abund),.packages = c('cluster','foreach'), .export = 'communityMST', .options.snow = opts) %dopar%
      communityMST(abund[i,], std_pool_traits, abundWeighted = abundWeighted, dis_metric = dis_metric)
    snow::stopCluster(cl)
  } else {
    comm_MST <- foreach(i=1:nrow(abund)) %do% {
      progress(i)
      communityMST(abund[i,], std_pool_traits, abundWeighted = abundWeighted, dis_metric = dis_metric)
    }
  }
  close(pb)
  
  message(sprintf("\nCalculating FEE0 for the %d communities...",nrow(abund)))
  pb <- utils::txtProgressBar(max = nrow(abund), initial = NA, style = 3)
  FEE0_comm <- function(mst_comm) {
    ll <- mst_comm
    mn_l <- mean(mst_comm)
    return(sum(sapply(ll,min,mn_l)))
  }
  if(doParallel & nrow(abund)>1) {
    cl <- snow::makeCluster(min(numCores,n_core_slurm,nrow(abund),na.rm=T))
    registerDoSNOW(cl)
    message(sprintf("With %d cores...",min(numCores,n_core_slurm,nrow(abund),na.rm=T)))
    progress(0)
    #registerDoParallel()
    opts <- list(progress = progress)
    FEE0 <- foreach(i=1:nrow(abund),.packages = c('cluster','foreach'), .options.snow = opts) %dopar% FEE0_comm(comm_MST[[i]])
    snow::stopCluster(cl)
  } else {
    FEE0 <- foreach(i=1:nrow(abund)) %do% {progress(i); FEE0_comm(comm_MST[[i]])}
  }
  FEE0 <- unlist(FEE0)
  close(pb)
  if(!is.null(user_ecdf)) {
    if(file.exists(user_ecdf)) {
      usr_d <- uni_d <- NULL
      load(user_ecdf)  # contains usr_d (pool-based eCDF), uni_d (uniform-based eCDF)
    } else message(sprintf("\nFile does not exist: %s",user_ecdf))
  }
  new_ecdf <- FALSE
  uni_nsp_comm <- unique(nsp_comm)
  if(1 %in% uni_nsp_comm) warning("\n  For monocultures, returning NA.")
  ### need to warn when a community has only 1 species
  
  if(poolBased) {
    message("\nPreparing pool-based eCDF...")
    if(!exists('usr_d',inherits=FALSE)) usr_d <- NULL
    fee0_ecdf <- foreach(n=uni_nsp_comm) %do% {   # cannot be %dopar": need to be serialized
      nm_ecdf <- sprintf('%s_%dd_%dsp',substr(dis_metric,1,3),num_traits,n)
      if(nm_ecdf%in%names(usr_d)) usr_d[[nm_ecdf]]        # use usr_d data if available
      else {
        new_ecdf <- TRUE    # new eCDF is generated
        if(n>1) message(sprintf("\nCumulative distribution function is unavailable for %d traits and %d species.\nIn calculating...",num_traits,n))
        n_ecdf <- FEE0_eCDF(n,num_traits,dis_metric=dis_metric,pool_traits=std_pool_traits)
        old_len <- length(usr_d)
        usr_d <- append(usr_d,list(n_ecdf))
        names(usr_d)[old_len+1] <- nm_ecdf
        n_ecdf
      }
    }
  } else {
    message("\nPreparing eCDF...")
    if(num_traits<=6) sys_d <- getFromNamespace(sprintf('fee_ecdf_%s_%dd',substr(dis_metric,1,3),num_traits),'FEE') # pre-calcuated eCDFs in sysdata.rda
    else sys_d <- NULL
    if(!exists('uni_d',inherits=FALSE)) uni_d <- NULL
    fee0_ecdf <- foreach(n=uni_nsp_comm) %do% {   # cannot be %dopar": need to be serialized
      nm_ecdf <- sprintf('%s_%dd_%dsp',substr(dis_metric,1,3),num_traits,n)
      if(nm_ecdf%in%names(sys_d)) sys_d[[nm_ecdf]]        # use sys_d data if available
      else if(nm_ecdf%in%names(uni_d)) uni_d[[nm_ecdf]]   # use uni_d data if available
      else {
        new_ecdf <- TRUE    # new eCDF is generated
        if(n>1) message(sprintf("\nCumulative distribution function is unavailable for %d traits and %d species.\nIn calculating...",num_traits,n))
        n_ecdf <- FEE0_eCDF(n,num_traits,dis_metric=dis_metric)
        old_len <- length(uni_d)
        uni_d <- append(uni_d,list(n_ecdf))
        names(uni_d)[old_len+1] <- nm_ecdf
        n_ecdf
      }
    }
  }
  if(new_ecdf&!is.null(user_ecdf)) {    # save usr_d and/or uni_d to file
    message(sprintf('\nSaving eCDF results to file: %s',user_ecdf))
    if(exists('usr_d',inherits=FALSE)&exists('uni_d',inherits=FALSE)) try(save(list=c('usr_d','uni_d'),file=user_ecdf))
    else if(exists('usr_d',inherits=FALSE)) try(save(list='usr_d',file=user_ecdf))
    else if(exists('uni_d',inherits=FALSE)) try(save(list='uni_d',file=user_ecdf))
  }
  names(fee0_ecdf) <- paste0('nsp_',uni_nsp_comm)
  message(sprintf("\nCalculating FEE for the %d communities...",nrow(abund)))
  pb <- utils::txtProgressBar(max = nrow(abund), initial = NA, style = 3)
  FEE <- foreach(i=1:nrow(abund)) %do% {
    progress(i)
    if(nsp_comm[i]<=1) NA
    else if(is.na(FEE0[i])) NA
    else fee0_ecdf[[paste0('nsp_',nsp_comm[i])]](FEE0[i])
  }
  close(pb)
  if(doFEE0) {
    return(list(FEE0=FEE0,FEE=unlist(FEE)))
  } else return(unlist(FEE))
}


#' Calculate species richness
#'
#' Computes species richness (taxonomic diversity) for the communities specified in \code{abund}.
#'
#' @param abund matrix (or data frame) of species abundances in the communities
#'   of interest. Columns are for different species in the species pool. Each
#'   row corresponds a community. For species not presented in a community,
#'   their abundance values should be assigned 0 or \code{NA}.
#' @return Number of species in each of the communities specified in \code{abund}.
#' @examples
#' d_comm <- rbind(comm1 = c(0:9),
#'                 comm2 = rep(1, 10),
#'                 comm3 = runif(10, 0, 5),
#'                 comm4 = rchisq(10, 2))
#' computeNSP(d_comm)
#' @export
computeNSP <- function(abund) {
  abund[is.na(abund)] <- 0
  return(apply(abund,1,function(x) length(which(x>0))))
}


# create the coordinates of n regular grids per each side in a unit hypercube with dimension of dim
coor_grids <- function(n,dim) {
  z <- data.frame(replicate(dim,seq(0,1,1/(n-1))))
  zz <- expand.grid(z)
  return(zz)
}

#' Calculate community's minimum spanning tree (MST)
#'
#' Get branch length information of MST for the community specified in
#' \code{abund} in the trait space defined by \code{spp_traits}.
#'
#' @param abund a numeric vector of species abundances in the community of
#'   interest. Species order in this vector needs to be same as that in species
#'   trait data (rows of \code{spp_traits}). Thus its length needs to be same as
#'   the number of rows in \code{spp_traits}.
#' @param spp_traits a matrix (or data frame) of numeric functional traits for
#'   species. Columns are for different traits and rows are for different
#'   species in the species pool.
#' @param abundWeighted logical; indicating whether the abundance values in
#'   \code{abund} are considered as weights in MST calculation. The default is
#'   \code{TRUE}.
#' @param dis_metric string specifying the scheme of quantifying species
#'   dissimilarity in trait space. The currently available options are
#'   \code{"euclidean"} (the default) and \code{"manhattan"}. Euclidean
#'   distances are root sum-of-squares of differences, and manhattan distances
#'   are the sum of absolute differences.
#' @details ...
#' @return If \code{abundWeighted} is \code{TRUE}, return sorted branch
#'   lengths of abundance-weighted MST. If \code{abundWeighted} is \code{FALSE},
#'   return sorted branch lengths of raw MST.
#' @examples
#' d_trait <- cbind(tr1 = rnorm(10, 0, 1), tr2 = rnorm(10, 3, 2))
#' comm1 = c(0:9)
#' comm2 = rep(1, 10)
#' communityMST(comm1, d_trait, abundWeighted = TRUE)
#' communityMST(comm1, d_trait)
#' communityMST(comm2, d_trait, abundWeighted = TRUE)
#' communityMST(comm2, d_trait)
#' @export
communityMST <- function(abund, spp_traits, abundWeighted = TRUE, dis_metric = c("euclidean","manhattan")) {
  dis_metric <- match.arg(dis_metric)
  abund <- as.vector(abund,mode="numeric")
  abund[is.na(abund)] <- 0
  if(length(abund)==0|sum(abund)==0) return(list(NA,NA))
  stopifnot(length(abund) == nrow(spp_traits))
  stopifnot(sum(abund) > 0)
  stopifnot(dis_metric %in% c("euclidean","manhattan"))
  x <- which(abund>0)
  w <- abund[x]
  nsp <- length(x)
  if(nsp<=1) return(list(NA,NA))
  m_dist <- dist(spp_traits[x,], method = dis_metric)
  if(length(unique(w))==1|!abundWeighted) {
    mst0 <- sort(m_dist[which(stats::as.dist(ape::mst(m_dist))==1)],decreasing=T)
    return(mst0)
  } else {
    p_w <- w/sum(w)
    m_w <- matrix(NA,nrow=nsp,ncol=nsp)
    rc <- combn(1:nsp,2)
    w <- combn(p_w,2,function(x) {
      k <- max(x)/min(x)
      min(1,k/(k+1)*2/nsp/mean(x))
    })
    m_w[t(rc)] <- w
    m_w[t(rc)[,2:1]] <- w
    m_dist_adj <- m_dist*as.dist(m_w)
    mst_adj <- sort(m_dist_adj[which(as.dist(ape::mst(m_dist_adj))==1)],decreasing=T)
    return(mst_adj)
  }
}


#' Create empirical cumulative distribution function of FEE0
#'
#' Use null model approach to generate empirical cumulative distribution function
#' (eCDF) of FEE0.
#'
#' @param num_species integer; number of species.
#' @param num_traits integer; number of traits (i.e., dimension of the trait
#' space).
#' @param dis_metric string specifying the scheme of quantifying species
#'   dissimilarity in trait space. The currently available options are
#'   \code{"euclidean"} (the default) and \code{"manhattan"}. Euclidean
#'   distances are root sum-of-squares of differences, and manhattan distances
#'   are the sum of absolute differences.
#' @param pool_traits a matrix (or data frame) of numeric functional traits for
#'   species pool. Columns are for different traits and rows are for different
#'   species in the species pool. If \code{pool_traits} is \code{NULL}, the
#'   calculation is then under the assumption of no prior knowledge of trait
#'   values (i.e., uniform distribution of trait values).
#' @return The eCDF of FEE0, used to transfer FEE0 metric to FEE index
#' @examples
#' d_trait <- cbind(tr1 = rnorm(10, 0, 1), tr2 = rnorm(10, 3, 2))
#' ecdf1 <- FEE0_eCDF(3, 4)
#' ecdf1(1.5)
#' ecdf2 <- FEE0_eCDF(9, 2, pool_traits = d_trait)
#' ecdf3 <- FEE0_eCDF(9, 2)
#' ecdf2(1.5)
#' ecdf3(1.5)
#' @export
FEE0_eCDF <- function(num_species, num_traits, dis_metric=c("euclidean","manhattan"), pool_traits=NULL) {
  if(num_species<=1) return(NA)
  message(sprintf("\nNumber of species: %d",num_species))
  progress <- function(n) utils::setTxtProgressBar(pb,n)
  if(is.null(pool_traits)) { # no prior knowledge of trait values (uniform distribution of trait values)
    pb <- utils::txtProgressBar(max = 10000, initial = NA, style = 3)
    nullMST <- foreach(i=1:10000) %do% {
      progress(i)
      coors <- replicate(num_traits,runif(num_species))
      if(num_traits==1) coors <- matrix(coors,ncol=1)
      communityMST(rep(1,num_species),coors,abundWeighted=FALSE,dis_metric=dis_metric)
    }
    close(pb)
  } else {  # pool-based eCDF
    if(num_traits!=ncol(pool_traits)) {
      warning('Numbers of traits are inconsistent...\n',
              sprintf('(num_trait: %d vs. columns of pool_traits: %d).\n',num_trait,ncol(pool_traits)),
              sprintf('Use column of pool_traits as the number of traits: %d',ncol(pool_traits)))
    }
    num_traits <- ncol(pool_traits)  # replace num_traits if in conflict
    nsp_pool <- nrow(pool_traits)
    if(!all(apply(pool_traits,2,max)==1)|!all(apply(pool_traits,2,min)==0)) pool_traits <- apply(pool_traits,2,function(x) (x-min(x))/(max(x)-min(x)))  # standardize traits to 0-1
    if(choose(nsp_pool,num_species)<=10000) {
      sp_ind <- combn(1:nsp_pool,num_species)
    }
    else sp_ind <- foreach(i=1:10000,.combine=cbind) %do% sample(nrow(pool_traits),num_species)
    pb <- utils::txtProgressBar(max = min(10000,choose(nsp_pool,num_species)), initial = NA, style = 3)
    nullMST <- foreach(i=1:min(10000,choose(nsp_pool,num_species))) %do% {
      progress(i)
      coors <- pool_traits[sp_ind[,i],]
      if(num_traits==1) coors <- matrix(coors,ncol=1)
      communityMST(rep(1,num_species),coors,abundWeighted=FALSE,dis_metric=dis_metric)
    }
    close(pb)
  }
  fee0 <- sapply(nullMST,function(x) sum(sapply(x,min,mean(x))))
  return(ecdf(fee0))
}
