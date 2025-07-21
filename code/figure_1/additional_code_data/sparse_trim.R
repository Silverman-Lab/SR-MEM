library(SparseDOSSA2)

generate_a <- function(n,
                       feature_param, Omega,
                       maxit = 10, verbose = FALSE) {
  i_iter <- 0
  samples_a <- matrix(NA,
                      nrow = nrow(feature_param),
                      ncol = 0)
  rownames(samples_a) <- rownames(feature_param)

  while(TRUE) {
    if(i_iter + 1 > maxit)
      stop("Can't satisfy conditions!")
    i_iter <- i_iter + 1
    if(verbose)
      print(i_iter)

    samples_a <- cbind(samples_a,
                       t(rcopulasso(n = n,
                                    pi0 = feature_param[, "pi0"],
                                    mu = feature_param[, "mu"],
                                    sigma = feature_param[, "sigma"],
                                    Omega = Omega)))

    ind_nonzero <- apply(samples_a > 0, 2, any) ## FIXME??
    if(sum(ind_nonzero) >= n) {
      samples_a <- samples_a[, ind_nonzero][, seq_len(n)]
      colnames(samples_a) <- paste0("Sample", seq_len(n))
      return(samples_a)
    }
  }
}

generate_count <- function(rel, depth) {
  mat_count <-
    vapply(seq_along(depth),
           function(i_sample) {
             if(all(rel[, i_sample] == 0))
               return(rep(0, nrow(rel)))
             else
               rmultinom(n = 1, size = depth[i_sample], prob = rel[, i_sample])
           },
           rep(0, nrow(rel)))
  dimnames(mat_count) <- dimnames(rel)
  return(mat_count)
}

rZILogN_one <- function(n, pi0, mu, sigma) {
    return(exp(rnorm(n = n, mean = mu, sd = sigma)) *
             rbinom(n = n, size = 1, prob = 1 - pi0))
}

TSS <- function(x, correct = FALSE) {
  if(any(x < 0))
    stop("Negative x values are not accepted!")
  if(all(x == 0)) return(x)
  # this is a special case where only one feature is present
  if(sum(x > 0) == 1 & correct) {
    x[x > 0] <- 1 - 0.5 / length(x)
    return(x)
  }
  return(x / sum(x))
}

generate_depth <- function(mu_depth,
                           sigma_depth,
                           n, median_depth) {
  depth <- exp(rnorm(n = n, mean = mu_depth, sd = sigma_depth))
  depth <- round(depth / median(depth) * median_depth)
  return(depth)
}

rcopulasso <- function(n, pi0, mu, sigma, Omega) {
  if(length(pi0) != length(mu) |
     length(pi0) != length(sigma) |
     length(pi0) != nrow(Omega))
    stop("Parameter dimensions must agree!")

  # sample marginals
  mat_amarginals <-
    vapply(seq_len(length(pi0)),
           function(i)
             rZILogN_one(n = n,
                         pi0 = pi0[i],
                         mu = mu[i],
                         sigma = sigma[i]),
           rep(0.0, n))
  # arrange from smallest to largest for shuffling
  mat_amarginals <-
    apply(mat_amarginals, 2, function(x) x[order(x)])

  # sample ranks
  mat_rank <-
    mvtnorm::rmvnorm(n = n, sigma = solve(Omega))
  mat_rank <- apply(mat_rank, 2, rank)

  mat_a <-
    vapply(seq_len(length(pi0)),
           function(i)
             mat_amarginals[, i, drop = TRUE][mat_rank[, i, drop = TRUE]],
           rep(0.0, n))

  return(mat_a)
}

logit <- function(x) log(x) - log(1 - x)

expit <- function(x)
  exp(x) / (1 + exp(x))

spike_oneA_metadata <- function(param,
                                metadata,
                                col_abundance = c(),
                                effect_abundance = c(),
                                col_prevalence = c(),
                                effect_prevalence = c()) {
  effect_abundance_all <- rep(0, ncol(metadata))
  effect_prevalence_all <- rep(0, ncol(metadata))

  effect_abundance_all[col_abundance] <- effect_abundance
  effect_prevalence_all[col_prevalence] <- effect_prevalence

  pi0 <- expit(logit(param["pi0"]) - (metadata %*% effect_prevalence_all)[, 1])
  mu <- param["mu"] + (metadata %*% effect_abundance_all)[, 1]

  a <- rZILogN_one(n = nrow(metadata),
                   pi0 = pi0,
                   mu = mu,
                   sigma = param["sigma"])
  return(a)
}

spike_a_metadata <- function(null,
                             feature_param,
                             metadata,
                             spike_df) {
  if(ncol(null) != nrow(metadata))
    stop("Sample size of null abundance matrix and metadata do not agree!")
  if(nrow(spike_df) !=
     nrow(dplyr::distinct(spike_df,
                          metadata_datum,
                          feature_spiked,
                          associated_property)))
    stop("feature-metadata spiking specification data frame has duplicate ",
         "metadata_datum, feature_spiked, and associated_property tuples!")
  spiked <- null
  for(feature_i in unique(spike_df$feature_spiked)) {
    spike_df_i <- subset(spike_df, feature_spiked == feature_i)
    spike_df_i_abundance <- subset(spike_df_i,
                                   associated_property == "abundance")
    spike_df_i_prevalence <- subset(spike_df_i,
                                    associated_property == "prevalence")
    spiked[feature_i, ] <-
      spike_oneA_metadata(param = feature_param[feature_i, ],
                          metadata = metadata,
                          col_abundance = spike_df_i_abundance$metadata_datum,
                          effect_abundance = spike_df_i_abundance$effect_size,
                          col_prevalence = spike_df_i_prevalence$metadata_datum,
                          effect_prevalence = spike_df_i_prevalence$effect_size)
  }

  dimnames(spiked) <- dimnames(null)
  return(spiked)
}

SparseDOSSA2_m <- function(template = "Stool",
                         n_sample = 100,
                         new_features = TRUE,
                         n_feature = 100,
                         spike_metadata = "none",
                         metadata_effect_size = 1,
                         perc_feature_spiked_metadata = 0.05,
                         metadata_matrix = NULL,
                         median_read_depth = 50000,
                         verbose = TRUE) {
  if(is.character(template)) {
    if(!template %in% c("Stool", "Vaginal", "IBD"))
      stop("Pre-trained template must be one of \"Stool\", \"Vaginal\", or \"IBD\"!")
    template <- get(template)
  }

  # generate per-feature params
  feature_param_template <- cbind("pi0" = template$EM_fit$fit$pi0,
                                  "mu" = template$EM_fit$fit$mu,
                                  "sigma" = template$EM_fit$fit$sigma)
  if(!new_features) {
    if(verbose)
      message("new_features is FALSE, ",
              "adopt per-feature parameters from template ",
              "(n_feature will be ignored)...")
    n_feature <- nrow(feature_param_template)
    feature_param <- feature_param_template
    Omega <- template$EM_fit$fit$Omega
    # check that Omega and Sigma should agree
    if(max(abs(Omega %*% template$EM_fit$fit$Sigma -
               diag(rep(1, length(template$EM_fit$fit$pi0))))) > 1e-5)
       stop("Omega shoud be the inverse of Sigma!")
  } else {
    if(verbose)
      message("new_features is TRUE, ",
              "generating new per-feature parameters, based on template...")
    feature_param <-
      generate_featureParam(F_fit = template$F_fit,
                            param_original = feature_param_template,
                            n_feature = n_feature)
    Omega <- diag(rep(1, nrow(feature_param)))
  }
  features <-
    rownames(Omega) <-
    colnames(Omega) <-
    rownames(feature_param)

  # generate null absolute abundance matrix
  if(verbose)
    message("Generating null absolute abundance matrix...")
  mat_null <- generate_a(n = n_sample,
                         feature_param = feature_param,
                         Omega = Omega)

  # generate spiked-in association with metadata
  if(!(is.character(spike_metadata) | is.data.frame(spike_metadata)))
    stop("spike_metadata must be either character or a data.frame!")
  if(is.character(spike_metadata)) {
    spike_metadata <- match.arg(spike_metadata,
                                choices = c("none", "abundance",
                                            "prevalence", "both"))
  }
  if(identical(spike_metadata, "none")) {
    if(verbose)
      message("spike_metadata is \"none\", ",
              "no metadata association will be simulated...")

    mat_spiked_metadata <- mat_null
    feature_metadata_spike_df <- NULL
  } else {
    if(verbose)
      message("Spiking in metadata association...")

    # metadata_matrix
    if(is.null(metadata_matrix)) {
      if(is.data.frame(spike_metadata))
        stop("spike_metadata is provided as a data.frame. User must specify ",
             "metadata_matrix as well!")
      if(verbose)
        message("metadata_matrix is not provided; ",
                "simulating default metadata_matrix...")
      metadata_matrix <- cbind(rnorm(n = n_sample),
                               rbinom(n = n_sample,
                                      size = 1,
                                      prob = 0.5))
      rownames(metadata_matrix) <- colnames(mat_null)
    } else {
      if(!is.matrix(metadata_matrix))
        stop("metadata_matrix must be a matrix ",
             "(model matrix where categorical variables are dummified)!")
      if(nrow(metadata_matrix) != n_sample)
        stop("n_sample does not agree with number of samples in ",
             "metadata_matrix!")
      if(!is.null(rownames(metadata_matrix)))
        colnames(mat_null) <- rownames(metadata_matrix)
    }
    n_metadata <- ncol(metadata_matrix)

    # metadata_effect_size
    if(length(metadata_effect_size) != 1 &
       length(metadata_effect_size) != n_metadata)
      stop("Length of metadata_effect_size can only be either 1 or number of ",
           "columns of metadata_matrix!")
    if(length(metadata_effect_size) == 1)
      metadata_effect_size <- rep(metadata_effect_size, n_metadata)

    # feature_metadata_spike_df
    if(is.data.frame(spike_metadata)) {
      if(verbose) {
        message("spike_metadata is provided as a data.frame; ",
                "will use for simulating metadata association ",
                "(metadata_effect_size and perc_feature_spiked_metadata will ",
                "be ignored)...")
      }

      # check format
      if(!all(c("metadata_datum",
                "feature_spiked",
                "associated_property",
                "effect_size") %in%
              colnames(spike_metadata)))
        stop("spike_metadata does not follow the correct format! ",
             "Must have the following columns: metadata_datum, ",
             "feature_spiked, associated_property, and effect_size.")
      if(!all(spike_metadata$feature_spiked %in%
              features))
        stop("feature_spiked in spike_metadata must provide the ",
             "spiked feature names!")
      if(!all(spike_metadata$metadata_datum %in%
              seq(1, n_metadata)))
        stop("metadata_datum in spike_metadata must provide the ",
             "associated metadata column number!")
      if(!all(spike_metadata$associated_property %in%
              c("prevalence", "abundance")))
        stop("associated_property in spike_metadata must be ",
             "either \"prevalence\" or \"abundance\"!")

      feature_metadata_spike_df <- spike_metadata
    } else {
      if(verbose)
        message("spike_metadata is specified as ", spike_metadata, "; ",
                "generating default metadata association...")
      feature_metadata_spike_df <-
        generate_feature_metadata_spike_df(
          features = features,
          perc_feature_spiked_metadata = perc_feature_spiked_metadata,
          n_metadata = n_metadata,
          effect_size = metadata_effect_size,
          spike_metadata = spike_metadata)
    }
    if(verbose)
      message("Generating feature abundances with spiked-in metadata ",
              "associations...")
    mat_spiked_metadata <-
      spike_a_metadata(null = mat_null,
                       feature_param = feature_param,
                       metadata = metadata_matrix,
                       spike_df = feature_metadata_spike_df)
  }

  mat_rel <- apply(mat_spiked_metadata, 2, TSS)
  if(any(is.na(template$depth_fit))) {
    mat_count <- mat_rel
  } else {
    if(verbose)
      message("Generating count matrix...")
    # generate read depth
    depth_new <- generate_depth(mu_depth = template$depth_fit["mu_depth"],
                                sigma_depth = template$depth_fit["sigma_depth"],
                                n = n_sample,
                                median_depth = median_read_depth)
    # generate read counts
    mat_count <- generate_count(rel = mat_rel,
                                depth = depth_new)
  }

  return(list(simulated_data = mat_count,
              simulated_matrices = list(rel = mat_rel,
                                        a_spiked = mat_spiked_metadata,
                                        a_null = mat_null),
              params = list(feature_param = feature_param,
                            Omega = Omega),
              template = template,
              spike_metadata = list(spike_metadata = spike_metadata,
                                    metadata_matrix = metadata_matrix,
                                    feature_metadata_spike_df =
                                      feature_metadata_spike_df)))
}
