# Adapted from https://enigma.ini.usc.edu/wp-content/uploads/combat_for_ENIGMA_sMRI/combat_for_ENIGMA_sMRI.R
# https://doi.org/10.1016/j.neuroimage.2020.116956

################################################################################
#                                                                              #
# combat_fit and combat_apply functions                                        #
# =====================================                                        #
#                                                                              #
# Authors:                                                                     #
#                                                                              #
#  1) The original ComBat function was in the sva package that can be found at #
#     https://bioconductor.org/packages/release/bioc/html/sva.html             #
#                                                                              #
#  2) First modification by Jean-Philippe Fortin for the harmonization of MRI  #
#     data. Please cite https://10.1016/j.neuroimage.2017.08.047               #
#                                                                              #
#  3) Second modification by Joaquim Radua to separate functions for fitting   #
#     and applying the harmonization, allow missings and constant rows and mi- #
#     nor changes in the arguments of the functions to facilitate their use.   #
#     Please cite "Increased power by harmonizing structural MRI site diffe-   #
#     rences with the ComBat batch adjustment method in ENIGMA".                #
#                                                                              #
# The original and present code is under the Artistic License 2.0. If using    #
# this code, make sure you agree and accept this license.                      #
#                                                                              #
################################################################################


# Calculate some characteristics on the batches
.combat_tmp1 <- function (dat, batch, levels_batch, mod) {
  batchmod <- model.matrix(~ -1 + batch)
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels_batch[i])
  }
  # List of samples in each batch
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  # Combine batch variable and covariates
  design <- cbind(batchmod, mod)
  # Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function (x) all(x == 1))
  design <- as.matrix(design[, !check])
  batch.design <- design[,1:n.batch]
  return(list(
    dat = dat,
    batchmod = batchmod,
    n.batch = n.batch,
    batches = batches,
    n.batches = n.batches,
    n.array = n.array,
    design = design,
    batch.design = batch.design
  ))
}
# Estimate B.hat, grand.mean and var.pooled
.combat_tmp2 <- function (tmp1, verbose = TRUE) {
  # Number of covariates or covariate levels
  if (verbose) {
    cat(
      "[combat] Adjusting for",
      ncol(tmp1$design) - ncol(tmp1$batchmod),
      "covariate(s) or covariate level(s)\n"
    )
  }
  # Check if the design is confounded
  if (qr(tmp1$design)$rank < ncol(tmp1$design)) {
    if (ncol(tmp1$design) == (tmp1$n.batch + 1)) {
      stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if (ncol(tmp1$design) > (tmp1$n.batch + 1)) {
      if ((qr(tmp1$design[,-c(1:tmp1$n.batch)])$rank < ncol(tmp1$design[,-c(1:tmp1$n.batch)]))) {
        stop("The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.")
      } else {
        stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
  }
  # Standardize data across features
  B.hat <- solve(t(tmp1$design) %*% tmp1$design) %*% t(tmp1$design) %*% t(as.matrix(tmp1$dat))
  # Standarization Model
  grand.mean <- t(tmp1$n.batches / tmp1$n.array) %*% B.hat[1:tmp1$n.batch,]
  var.pooled <- ((tmp1$dat - t(tmp1$design %*% B.hat))^2) %*% rep(1 / tmp1$n.array, tmp1$n.array)
  return(list(
    B.hat = B.hat,
    grand.mean = grand.mean,
    var.pooled = var.pooled
  ))
}
# Standardize data
.combat_tmp3 <- function (dat, tmp1, tmp2, verbose = TRUE) {
  if (verbose) {
    cat("[combat] Standardizing data across features\n")
  }
  stand.mean <- t(tmp2$grand.mean) %*% t(rep(1, tmp1$n.array))
  if (!is.null(tmp1$design)) {
    tmp <- tmp1$design;tmp[,c(1:tmp1$n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% tmp2$B.hat)
  } 
  s.data <- (dat-stand.mean) / (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))
  return(list(
    stand.mean = stand.mean,
    s.data = s.data
  ))
}
# Fit L/S model 
.combat_tmp4 <- function (tmp1, tmp2, tmp3, eb = TRUE, verbose = TRUE) {
  # Get regression batch effect parameters
  if (eb) {
    if (verbose) {
      cat("[combat] Fitting L/S model and finding priors\n")
    }
  } else {
    if (verbose) {
      cat("[combat] Fitting L/S model\n")
    }
  }
  gamma.hat <- solve(t(tmp1$batch.design) %*% tmp1$batch.design) %*% t(tmp1$batch.design) %*% t(as.matrix(tmp3$s.data))
  delta.hat <- NULL
  for (i in tmp1$batches) {
    delta.hat <- rbind(delta.hat, apply(tmp3$s.data[,i], 1, var, na.rm = T))
  }
  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb) {
    # Find Priors
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, .aprior)
    b.prior <- apply(delta.hat, 1, .bprior)
    # Find EB batch adjustments
    if (verbose) {
      cat("[combat] Finding parametric adjustments\n")
    }
    for (i in 1:tmp1$n.batch) {
      temp <- .it.sol(tmp3$s.data[,tmp1$batches[[i]]], gamma.hat[i,], delta.hat[i,], gamma.bar[i], t2[i], a.prior[i], b.prior[i])
      gamma.star <- rbind(gamma.star, temp[1,])
      delta.star <- rbind(delta.star, temp[2,])
    }
  } 
  return(list(
    gamma.hat = gamma.hat,
    delta.hat = delta.hat, 
    gamma.star = gamma.star,
    delta.star = delta.star, 
    gamma.bar = gamma.bar,
    t2 = t2,
    a.prior = a.prior,
    b.prior = b.prior
  ))
}
# Adjust the data
.combat_tmp5 <- function (tmp1, tmp2, tmp3, tmp4, eb = TRUE, verbose = TRUE) {
  # Normalize the data
  if (verbose) {
    cat("[combat] Adjusting the data\n")
  }
  bayesdata <- tmp3$s.data
  j <- 1
  for (i in tmp1$batches) {
    if (eb) {
      bayesdata[,i] <- (bayesdata[,i] - t(tmp1$batch.design[i,] %*% tmp4$gamma.star)) / (sqrt(tmp4$delta.star[j,]) %*% t(rep(1, tmp1$n.batches[j])))
    } else {
      bayesdata[,i] <- (bayesdata[,i] - t(tmp1$batch.design[i,] %*% tmp4$gamma.hat )) / (sqrt(tmp4$delta.hat[j,] ) %*% t(rep(1, tmp1$n.batches[j])))
    }
    j <- j + 1
  }
  return((bayesdata * (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))) + tmp3$stand.mean)
}
# Fit the combat harmonization model
combat_fit <- function (dat, batch, mod = NULL, eb = TRUE, verbose = TRUE) {
  # Modified by Joaquim Radua
  if (!is.factor(batch)) {
    stop("batch must be a factor")
  }
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  } else if (!is.matrix(dat)) {
    stop("dat must be a matrix")
  }
  if (ncol(dat) == length(batch)) {
    transpose <- FALSE
    if (verbose) {
      cat("[combat] Subjects are COLUMNS\n")
    }
  } else if (nrow(dat) == length(batch)) {
    transpose <- TRUE
    dat <- t(dat)
    if (verbose) {
      cat("[combat] Subjects are ROWS\n")
    }
  } else {
    stop("dat must have the same number of columns or rows than the length of batch")
  }
  if (is.data.frame((mod))) {
    mod <- as.matrix(mod)
  } else if (!(is.matrix(mod) || is.null(mod))) {
    stop("mod must be a matrix or NULL")
  }
  # Imputation of missing values
  if (any(is.na(dat))) {
    if (verbose) {
      cat("[combat] Imputing missing data (only for fit)\n")
    }
    for (batch_i in sort(unique(batch))) {
      i <- which(batch == batch_i)
      dat_i <- dat[,i]
      if (any(is.na(dat_i))) {
        for (j in 1:nrow(dat)) {
          dat_ji <- dat_i[j,]
          is_na <- which(is.na(dat_ji))
          # Some missing, impute from other subjects of the site
          if (length(is_na) > 0 && length(is_na) < length(i)) {
            if (length(is_na) == 1) {
              mod_i_is_na <- matrix(mod[i[is_na],], nrow = 1)
            } else {
              mod_i_is_na <- mod[i[is_na],]
            }
            beta <- matrix(coef(lm(dat_ji ~ mod[i,])))
            beta[which(is.na(beta))] <- 0
            dat[j,i[is_na]] <- cbind(1, mod_i_is_na) %*% beta
          } else {
            dat[j,i[is_na]] <- mean(dat_ji, na.rm = TRUE)
          }
        }
      }
    }
    # If there are still missing data, they are because a site has all missing
    # data for a ROI
    for (batch_i in sort(unique(batch))) {
      i <- which(batch == batch_i)
      if (any(is.na(dat[,i]))) {
        for (j in 1:nrow(dat)) {
          dat_j <- dat[j,]
          if (is.na(dat_j[i[1]])) {
            if (!is.null(mod)) {
              beta <- matrix(coef(lm(dat_j ~ mod)))
              beta[which(is.na(beta))] <- 0
              dat[j,i] <- cbind(1, mod[i,]) %*% beta
            } else {
              dat[j,i] <- mean(dat_j, na.rm = TRUE)
            }
          }
        }
      }
    }
  }
  not_constant <- which(apply(dat, 1, function (x) {var(x) > 0}))
  dat <- dat[not_constant,]
  if (eb) {
    if (verbose) {
      cat("[combat] Performing ComBat with empirical Bayes\n")
    }
  } else {
    if (verbose) {
      cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
    }
  }
  if (verbose) {
    cat("[combat] Found", nlevels(batch), "batches (e.g., sites)\n")
  }
  levels_batch <- levels(batch)
  tmp1 <- .combat_tmp1(dat, batch, levels_batch, mod)
  tmp2 <- .combat_tmp2(tmp1, verbose)
  tmp3 <- .combat_tmp3(dat, tmp1, tmp2, verbose)
  tmp4 <- .combat_tmp4(tmp1, tmp2, tmp3, eb, verbose)
  return(list(
    levels_batch = levels_batch,
    transpose = transpose,
    not_constant = not_constant,
    eb = eb,
    tmp2 = tmp2,
    tmp4 = tmp4
  ))
}
# Apply a combat harmonization model. The first argument is the result of
# "combat_fit"
combat_apply <- function (tmp, dat, batch, mod = NULL, verbose = TRUE) {
  # Modified by Joaquim Radua
  if (!is.factor(batch) || nlevels(batch) != length(tmp$levels_batch) || levels(batch) != tmp$levels_batch) {
    stop("batch must be a factor with the same levels than when fitting combat")
  }
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  } else if (!is.matrix(dat)) {
    stop("dat must be a matrix")
  }
  if (tmp$transpose) {
    if (nrow(dat) != length(batch)) {
      stop("dat must have the same number of rows than the length of batch")
    }
    dat <- t(dat)
  } else {
    if (ncol(dat) != length(batch)) {
      stop("dat must have the same number of columns than the length of batch")
    }
  }
  dat.combat <- dat
  dat <- dat[tmp$not_constant,]
  if (is.data.frame((mod))) {
    mod <- as.matrix(mod)
  } else if (!(is.matrix(mod) || is.null(mod))) {
    stop("mod must be a matrix or NULL")
  }
  tmp1 <- .combat_tmp1(dat, batch, tmp$levels_batch, mod)
  tmp3 <- .combat_tmp3(dat, tmp1, tmp$tmp2, verbose)
  tmp5 <- .combat_tmp5(tmp1, tmp$tmp2, tmp3, tmp$tmp4, tmp$eb, verbose)
  dat.combat[tmp$not_constant,] <- tmp5
  if (tmp$transpose) {
    dat.combat <- t(dat.combat)
  }
  return(list(
    dat.combat = dat.combat, 
    gamma.hat = tmp$tmp4$gamma.hat,
    delta.hat = tmp$tmp4$delta.hat, 
    gamma.star = tmp$tmp4$gamma.star,
    delta.star = tmp$tmp4$delta.star, 
    gamma.bar = tmp$tmp4$gamma.bar,
    t2 = tmp$tmp4$t2,
    a.prior = tmp$tmp4$a.prior,
    b.prior = tmp$tmp4$b.prior,
    batch = tmp1$batch,
    mod = tmp1$mod, 
    stand.mean = tmp3$stand.mean,
    stand.sd = sqrt(tmp$tmp2$var.pooled)[,1]
  ))
}


################################################################################
#                                                                              #
# This is a copy of the original code from the standard version of the sva     #
# package that can be found at                                                 #
# https://bioconductor.org/packages/release/bioc/html/sva.html                 #
# The original and present code is under the Artistic License 2.0.             #
# If using this code, make sure you agree and accept this license.             #
#                                                                              #
################################################################################


# Following four find empirical hyper-prior values
.aprior <- function (gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2 * s2 + m^2) / s2
}
.bprior <- function (gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m * s2 + m^3) / s2
}
.postmean <- function (g.hat, g.bar, n, d.star, t2) {
  (t2 * n * g.hat + d.star * g.bar) / (t2 * n + d.star)
}
.postvar <- function (sum2, n, a, b) {
  (0.5 * sum2 + b) / (n / 2 + a - 1)
}
# Pass in entire data set, the design matrix for the entire data, the batch
# means, the batch variances, priors (m, t2, a, b), columns of the data 
# matrix for the batch. Uses the EM to find the parametric batch adjustments
.it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) {
  n <- apply(!is.na(sdat), 1, sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- .postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 1, sum, na.rm = T)
    d.new <- .postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  # cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}
