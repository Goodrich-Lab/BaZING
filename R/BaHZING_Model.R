#' BaHZING_Model Function
#' This function implements the BaHZING model for microbiome data analysis.
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr
#' @importFrom utils globalVariables
#' @importFrom stats quantile update
#' @importFrom bayestestR p_direction p_rope p_map
#' @importFrom glue glue
#' @param formatted_data An object containing formatted microbiome data.
#' @param x A vector of column names of the exposures.
#' @param covar An optional vector of the column names of covariates.
#' @param n.chains An optional integer specifying the number of parallel chains
#' for the model in jags.model function. Default is 3.
#' @param n.adapt An optional integer specifying the number of iterations for
#' adaptation in jags.model function. Default is 100.
#' @param n.iter.burnin An optional integer specifying number of iterations in
#' update function. Default is 1,000.
#' @param n.iter.sample An optional integer specifying the number of iterations
#' in coda.samples function. Default is 5,000.
#' @param exposure_standardization Method for standardizing the exposures.
#' Should be one of "standard_normal" (the default), "quantile", or "none". If
#' "none", exposures are not standardized before analysis, and counterfactual
#' profiles must be specified by the user.
#' @param counterfactual_profiles A 2xP matrix or a vector with length of 2; P
#' is the number of exposures in x. If a 2xP matrix is provided,
#' the effect estimates for the mixture are interpreted as the estimated change
#' in the outcome when changing each exposure p in 1:P is changed from
#' `counterfactual_profiles[1,p]` to `counterfactual_profiles[2,p]`. If a vector of
#' length 2 is provided, the effect estimates for the mixture are interpreted as
#' the estimated change in the outcome when changing each exposure from
#' `counterfactual_profiles[1]` to `counterfactual_profiles[2]`. If
#' exposure_standardization = "standard_normal", then the default is c(-0.5, 0.5),
#' and the effect estimate is calculated based on increasing all exposures in
#' the mixture by one standard deviation. If exposure_standardization = "quantile",
#' then the default is c(0,1), and the effect estimate is calculated based on
#' increasing all exposures in the mixture by one quantile (where the number of
#' quantiles is based on the parameter q).
#' @param q An integer specifying the number of quantiles. Only required if
#' exposure_standardization = "quantile". If exposure_standardization =
#' "quantile" and q is not specified, then a default of q = 4 is used.
#' @param verbose If TRUE (default), function returns information a data quality
#' check.
#' @param return_all_estimates If FALSE (default), results do not include
#' the dispersion and omega estimates from the BaHZING model.
#' @param seed Optional positive integer used to generate reproducible,
#' distinct JAGS RNG initial values for each chain.
#' @param parallel Logical; if TRUE, chains can be run in parallel with
#' `parallel::mclapply()` on supported non-Windows platforms. Windows uses
#' sequential chains.
#' @param n.cores Number of cores to use when parallel is TRUE. If NULL,
#' physical cores are detected automatically. The value is capped at
#' `n.chains`.
#' @param ROPE_range Region of practical equivalence (ROPE) for calculating
#' p_rope. Default is c(-0.1, 0.1).
#' @param return_posterior If TRUE, the full posterior MCMC samples are
#' returned alongside the results data frame. Default is FALSE.
#' @param offset Logical; if TRUE (default), include the log of each sample's
#' total library size as an offset in both the count and zero-inflation model
#' predictors. If FALSE, no library-size offset is used.
#' @return A data frame containing results of the Bayesian analysis, with the
#' following columns:
#' - taxa_full: Full Taxa information, including all levels of the taxonomy.
#' Taxanomic levels are split by two underscores ('__').
#' - taxa_name: Taxa name, which is the last level of the taxonomy.
#' - domain: domain of the taxa.
#' - exposure: Exposure name (either one of  the individual exposures, or the
#' mixture).
#' - component: Zero inflated model estimate or the Count model estimate.
#' - estimate: Point estimate of the posterior distributions.
#' - bci_lcl: 95% Bayesian Credible Interval Lower Limit. Calculated as the
#' equal tailed interval of posterior distributions using the quantiles method.
#' - bci_ucl: 95% Bayesian Credible Interval Upper Limit. Calculated as the
#' equal tailed interval of posterior distributions using the quantiles method.
#' - p_direction: The Probability of Direction, calculated with `bayestestR`. A
#' higher value suggests a higher probability that the estimate is strictly
#' positive or negative. In other words, the closer the value to 1, the higher
#' the probability that the estimate is non-zero. Values can not be less than
#' 50%. From `bayestestR`: also known as the Maximum Probability of Effect
#' (MPE). This can be interpreted as the probability that a parameter (described
#' by its posterior distribution) is strictly positive or negative (whichever
#' is the most probable). Although differently expressed, this index is fairly
#' similar (i.e., is strongly correlated) to the frequentist p-value.
#' - p_rope: The probability that the estimate is not within the Region of
#' practical equivalence (ROPE), calculated with `bayestestR`. The proportion
#' of the whole posterior distribution that doesn't lie within the `ROPE_range`.
#' - p_map: Bayesian equivalent of the p-value, calculated with `bayestestR`.
#'  From `bayestestR`:  p_map is related to the odds that a parameter (described
#'  by its posterior distribution) has against the null hypothesis (h0) using
#'   Mills' (2014, 2017) Objective Bayesian Hypothesis Testing framework. It
#'   corresponds to the density value at the null (e.g., 0) divided by the
#'   density at the Maximum A Posteriori (MAP).
#' @export
#' @name BaHZING_Model

# Declare global variables
globalVariables(c("X2.5.", "X97.5.", "Mean",
                  "Exposure.Index", "taxa_index", "taxa_full",
                  "component", "estimate", "bci_lcl", "bci_ucl",
                  "domain", "taxa_name", "pdir","prope","pmap",
                  "name"))

# JAGS model text generators -------------------------------------------------
# BaHZING's hierarchical shrinkage model has one recurring pattern: every
# taxonomic level's exposure effect is drawn from a normal distribution
# centered on its parent level's corresponding effect, borrowed via a binary
# taxonomy indicator matrix (inprod against e.g. GenusData). The narrowest
# level in taxa_levels (species, by default) is the one exception with a
# parent to borrow from AND its own data likelihood; the broadest level
# (phylum, by default) is the other exception, with no parent to borrow from
# at all. These generators build the JAGS text for an arbitrary-length
# taxa_levels hierarchy from these three roles (narrowest/data-likelihood,
# mid-hierarchy, broadest/terminal), replacing what used to be two ~180-line
# hand-maintained, near-identical literal JAGS strings (one per
# with/without-covariates variant) that assumed exactly 6 fixed level names.
#
# All generators use glue::glue(.open="<<", .close=">>", .trim=FALSE): custom
# delimiters because JAGS's own "{ }" block syntax would otherwise collide
# with glue's default interpolation markers, and .trim=FALSE because glue's
# default (.trim=TRUE) strips leading whitespace and the trailing newline -
# harmless for a block that ends the whole template, but it silently glues
# adjacent fragments onto the same line wherever one block's output is
# spliced into another.

#' Derive the standard JAGS identifier names for one taxa_levels position
#'
#' Single source of truth for the level-name -> JAGS-identifier convention
#' used throughout this file (mirrors `hierarchy_matrix_name()`'s role for
#' the incidence-matrix naming convention in `taxa_levels.R`): a level's
#' lowercase variable-name prefix, its proper-case label (for comments/
#' taxon-count scalar names), its loop index variable, its taxon-count
#' scalar name, and its incidence-matrix data variable name.
#'
#' @param taxa_levels Character vector, broadest to narrowest.
#' @param i Integer position to look up.
#' @return A list with elements `level` (lowercase prefix, e.g. `"genus"`),
#'   `label` (proper case, e.g. `"Genus"`), `idx` (loop index variable, e.g.
#'   `"g.r"`), `R` (taxon-count scalar name, e.g. `"Genus.R"`), and `data`
#'   (incidence-matrix variable name, e.g. `"GenusData"`).
#' @keywords internal
#' @noRd
.bahzing_level_vars <- function(taxa_levels, i) {
  label <- taxa_level_name(taxa_levels, i)
  level <- tolower(label)
  list(
    level = level,
    label = label,
    idx   = paste0(substr(level, 1, 1), ".r"),
    R     = paste0(label, ".R"),
    data  = paste0(label, "Data")
  )
}

#' Build the precision-prior and g-estimation lines for one taxonomic level
#'
#' Generates the JAGS lines for one level's precision priors (count-model and
#' zero-inflation dispersion) and its g-estimation block (the mixture
#' contrast computed from the low/high counterfactual `profiles`). Shared by
#' every level via `.bahzing_level_block()` and `.bahzing_narrowest_level_block()`.
#'
#' @param level Character. The level's variable-name prefix, e.g. `"genus"`.
#'   Used for the g-estimation nodes (`<level>.psi`, `<level>.eta.*`), which
#'   are always named after the level itself.
#' @param idx Character. The JAGS loop index variable for this level, e.g.
#'   `"g.r"`.
#' @return A length-1 character string ending in a trailing newline.
#' @keywords internal
#' @noRd
.bahzing_precision_gestimation_block <- function(level, idx) {
  tau_prefix <- paste0(level, ".")
  glue::glue(
"      # prior on precision
      <<tau_prefix>>tau[<<idx>>] <- 1/(<<tau_prefix>>sigma[<<idx>>]*<<tau_prefix>>sigma[<<idx>>])
      <<tau_prefix>>sigma[<<idx>>] ~ dunif(0,3)
      <<tau_prefix>>tau.zero[<<idx>>] <- 1/(<<tau_prefix>>sigma.zero[<<idx>>]*<<tau_prefix>>sigma.zero[<<idx>>])
      <<tau_prefix>>sigma.zero[<<idx>>] ~ dunif(0,3)

      # g-estimation
      <<level>>.eta.low[<<idx>>] <- inprod(<<level>>.beta[<<idx>>,1:P], profiles[1,1:P])
      <<level>>.eta.high[<<idx>>] <- inprod(<<level>>.beta[<<idx>>,1:P], profiles[2,1:P])
      <<level>>.psi[<<idx>>] <- <<level>>.eta.high[<<idx>>]-<<level>>.eta.low[<<idx>>]
      # zero-inflation
      <<level>>.eta.low.zero[<<idx>>] <- inprod(<<level>>.beta.zero[<<idx>>,1:P], profiles[1,1:P])
      <<level>>.eta.high.zero[<<idx>>] <- inprod(<<level>>.beta.zero[<<idx>>,1:P], profiles[2,1:P])
      <<level>>.psi.zero[<<idx>>] <- <<level>>.eta.high.zero[<<idx>>]-<<level>>.eta.low.zero[<<idx>>]
", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build the exposure-effect prior for one taxonomic level
#'
#' Generates `for(p in 1:P) { <level>.beta[idx,p] ~ dnorm(...) }`: each
#' level's exposure effect is drawn from a normal distribution centered on a
#' value borrowed from the parent level, computed via `inprod()` over the
#' parent's beta and the taxonomy indicator matrix - both the count-model and
#' zero-inflation components. If `parent` is `NULL` (only the broadest/
#' terminal level, phylum by default), there's nothing to borrow from, so the
#' prior falls back to a fixed `dnorm(0, tau)` instead.
#'
#' Shared by every level including the narrowest level, via
#' `.bahzing_level_block()` and `.bahzing_narrowest_level_block()`
#' respectively - both the borrowed-mean node name (`mu.<parent>`) and the
#' precision-node prefix (`<level>.tau`/`<level>.sigma`) now follow the same
#' convention for every level with no per-role override, unlike BaHZING's
#' original hand-written model (which named the narrowest level's precision
#' nodes bare `tau`/`sigma` and its borrowed-mean node `mu.species` instead
#' of `mu.genus`) - safe to make consistent since none of `tau`/`sigma`/
#' `mu.*` are ever monitored/extracted nodes.
#'
#' @param level Character. This level's variable-name prefix (e.g.
#'   `"genus"`), used for `<level>.beta`/`<level>.tau`.
#' @param idx Character. This level's JAGS loop index variable (e.g.
#'   `"g.r"`).
#' @param parent Character or `NULL`. The parent level's variable-name prefix
#'   to borrow the prior mean from (e.g. `"family"` for genus). `NULL` only
#'   for the broadest/terminal level, which has no parent.
#' @param parent_R Character or `NULL`. The JAGS scalar holding the parent
#'   level's taxon count (e.g. `"Family.R"`), used to size the `inprod()`.
#'   Required whenever `parent` is supplied.
#' @param parent_data Character or `NULL`. The taxonomy indicator matrix
#'   linking this level's taxa to the parent's (e.g. `"FamilyData"`).
#'   Required whenever `parent` is supplied.
#' @return A length-1 character string ending in a trailing newline.
#' @keywords internal
#' @noRd
.bahzing_exposure_prior_block <- function(level, idx, parent = NULL, parent_R = NULL,
                                            parent_data = NULL) {
  tau_prefix <- paste0(level, ".")
  body <- if (!is.null(parent)) {
    glue::glue(
"        <<level>>.beta[<<idx>>,p] ~ dnorm(mu.<<parent>>[<<idx>>,p], <<tau_prefix>>tau[<<idx>>])
        mu.<<parent>>[<<idx>>,p] <- inprod(<<parent>>.beta[1:<<parent_R>>,p], <<parent_data>>[<<idx>>,1:<<parent_R>>])
        #Zero inflation component
        <<level>>.beta.zero[<<idx>>,p] ~ dnorm(mu.<<parent>>.zero[<<idx>>,p], <<tau_prefix>>tau.zero[<<idx>>])
        mu.<<parent>>.zero[<<idx>>,p] <- inprod(<<parent>>.beta.zero[1:<<parent_R>>,p], <<parent_data>>[<<idx>>,1:<<parent_R>>])
", .open = "<<", .close = ">>", .trim = FALSE)
  } else {
    glue::glue(
"        <<level>>.beta[<<idx>>,p] ~ dnorm(0, <<tau_prefix>>tau[<<idx>>])
        #Zero inflation component
        <<level>>.beta.zero[<<idx>>,p] ~ dnorm(0, <<tau_prefix>>tau.zero[<<idx>>])
", .open = "<<", .close = ">>", .trim = FALSE)
  }
  glue::glue(
"      for(p in 1:P) {
<<body>>      }
", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build one mid-hierarchy or broadest/terminal level's full JAGS block
#'
#' Assembles one level's complete `for(idx in 1:level_R) { ... }` block: a
#' comment header, the exposure-effect prior
#' (`.bahzing_exposure_prior_block()`), and the precision-prior/g-estimation
#' lines (`.bahzing_precision_gestimation_block()`). Handles both mid-
#' hierarchy levels (`i > 1` - borrows its prior mean from the parent level
#' at `i - 1`) and the broadest/terminal level (`i == 1` - the top of the
#' hierarchy, with nothing to borrow from).
#'
#' This level's and its parent's JAGS identifiers (variable-name prefix,
#' loop index, taxon-count scalar, incidence-matrix name) are derived
#' automatically from `taxa_levels` and `i` via `.bahzing_level_vars()` -
#' callers only need to supply the hierarchy and a position.
#'
#' Not used for the narrowest level - that level has its own
#' likelihood/dispersion/intercept/covariate structure with no equivalent at
#' any other level, so it gets its own function,
#' `.bahzing_narrowest_level_block()`, which reuses the same two shared
#' sub-block generators this function uses.
#'
#' @param taxa_levels Character vector naming the taxonomic hierarchy,
#'   ordered broadest to narrowest.
#' @param i Integer position of this level in `taxa_levels`.
#' @return A length-1 character string: this level's complete JAGS block.
#' @keywords internal
#' @noRd
.bahzing_level_block <- function(taxa_levels, i) {
  cur    <- .bahzing_level_vars(taxa_levels, i)
  parent <- if (i > 1) .bahzing_level_vars(taxa_levels, i - 1) else NULL

  exposure_prior   <- .bahzing_exposure_prior_block(cur$level, cur$idx,
                        parent$level, parent$R, parent$data)
  precision_gestim <- .bahzing_precision_gestimation_block(cur$level, cur$idx)

  glue::glue(
"    # <<cur$label>> level
    for(<<cur$idx>> in 1:<<cur$R>>) {
<<exposure_prior>><<precision_gestim>>    }

", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build the narrowest level's JAGS block - the data-likelihood level
#'
#' The narrowest level in `taxa_levels` is tied directly to the observed
#' data, so unlike `.bahzing_level_block()` (mid-hierarchy/terminal levels,
#' which are pure prior/hierarchy blocks), this one also includes: the
#' zero-inflated negative binomial likelihood (`Y ~ dnegbin(...)`,
#' `zero ~ dbern(...)`), the dispersion prior (`disp`), the intercept priors
#' (`alpha`/`alpha.zero`), and - when `has_covar` is `TRUE` - the covariate
#' terms and their priors (`delta`/`delta.zero`). This is what collapses the
#' original model's separate with/without-covariates literal text into one
#' function.
#'
#' As with `.bahzing_level_block()`, this level's and its parent's JAGS
#' identifiers are derived automatically from `taxa_levels` and `i` via
#' `.bahzing_level_vars()` - the parent is always `i - 1` (`n >= 2` is
#' enforced by `validate_taxa_levels()`, so a parent always exists here).
#'
#' Reuses the same two shared sub-block generators every other level uses,
#' called exactly the same way as `.bahzing_level_block()` calls them - this
#' level's borrowed-mean and precision-node naming now follows the same
#' convention as every other level (see `.bahzing_exposure_prior_block()`'s
#' docs for why that's safe, a deliberate departure from BaHZING's original
#' hand-written model's species-specific naming).
#'
#' @param taxa_levels Character vector naming the taxonomic hierarchy,
#'   ordered broadest to narrowest.
#' @param i Integer position of the narrowest level in `taxa_levels` (i.e.
#'   `length(taxa_levels)`).
#' @param has_covar Logical. Whether covariates are included in the model.
#'   When `TRUE`, adds the `delta`/`delta.zero` covariate terms to the
#'   likelihood and their `dnorm(0, 1.0E-02)` priors; when `FALSE`, omits
#'   them entirely (matching the original model's with/without-covariates
#'   variants). Confirmed by diffing the two original literal model strings:
#'   covariate terms only ever appear in this block - every other level's
#'   text is identical between the with/without-covariates variants.
#' @param has_offset Logical. Whether to add `log(L[i,1])` to the count and
#'   zero-inflation predictors.
#' @return A length-1 character string: the complete data-likelihood-level
#'   JAGS block.
#' @keywords internal
#' @noRd
.bahzing_narrowest_level_block <- function(taxa_levels, i, has_covar,
                                            has_offset = TRUE) {
  cur    <- .bahzing_level_vars(taxa_levels, i)
  parent <- .bahzing_level_vars(taxa_levels, i - 1)
  level  <- cur$level

  covar_lambda <- if (has_covar) " + inprod(delta[r, 1:Q], W[i,1:Q])" else ""
  covar_pi     <- if (has_covar) " + inprod(delta.zero[r, 1:Q], W[i,1:Q])" else ""
  offset_term  <- if (has_offset) " + log(L[i,1])" else ""
  covar_prior  <- if (has_covar)
"
      # prior on covariate effects
      for(q in 1:Q) {
        delta[r,q] ~ dnorm(0, 1.0E-02)
        delta.zero[r,q] ~ dnorm(0, 1.0E-02)
      }
" else ""

  exposure_prior <- .bahzing_exposure_prior_block(level, "r", parent = parent$level,
                                                    parent_R = parent$R, parent_data = parent$data)
  precision_gestim <- .bahzing_precision_gestimation_block(level, "r")

  glue::glue(
"    for(r in 1:R) {
      for(i in 1:N) {
        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
        log(lambda[i,r]) <- alpha[r] + inprod(<<level>>.beta[r,1:P], X.q[i,1:P])<<covar_lambda>><<offset_term>>

        # zero-inflation
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(<<level>>.beta.zero[r,1:P], X.q[i,1:P])<<covar_pi>><<offset_term>>
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02)
      alpha.zero[r] ~ dnorm(0, 1.0E-02)
<<covar_prior>>
<<exposure_prior>>
<<precision_gestim>>    }

", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Assemble the complete JAGS model text for BaHZING_Model()
#'
#' The top-level entry point for the templating system: builds every
#' taxonomic level's block (the narrowest level via
#' `.bahzing_narrowest_level_block()`; every other level via
#' `.bahzing_level_block()`, each simply passed `taxa_levels` and its own
#' position - parent resolution happens automatically inside those
#' functions) and concatenates them narrowest-to-broadest into one complete
#' `model { ... }` string, matching the original hand-written text's
#' declaration order, ready to pass to `jags.model(file = textConnection(...))`.
#'
#' @param taxa_levels Character vector naming the taxonomic hierarchy,
#'   ordered broadest to narrowest (e.g. `default_taxa_levels`). Must have
#'   already been validated with `validate_taxa_levels()`.
#' @param has_covar Logical. Whether the model includes covariates - passed
#'   straight through to `.bahzing_narrowest_level_block()`, the only level
#'   that varies its text based on this (see that function for why).
#' @param has_offset Logical. Whether the narrowest-level likelihood includes
#'   the library-size offset.
#' @return A length-1 character string: the complete JAGS model text, from
#'   `model {` through the final closing `}`.
#' @keywords internal
#' @noRd
.bahzing_build_model_text <- function(taxa_levels, has_covar,
                                      has_offset = TRUE) {
  n <- length(taxa_levels)
  blocks <- vector("list", n)
  for (i in seq_len(n)) {
    blocks[[i]] <- if (i == n) {
      # Narrowest level: the data-likelihood role. Always has a parent
      # (n >= 2 is enforced by validate_taxa_levels()).
      .bahzing_narrowest_level_block(taxa_levels, i, has_covar, has_offset)
    } else {
      .bahzing_level_block(taxa_levels, i)
    }
  }
  # Emit narrowest -> broadest (reverse of taxa_levels' broadest-first
  # storage order) to match the original text's declaration order. JAGS
  # registers stochastic nodes in declaration order, which can affect the
  # sampler's RNG draw sequence for a given seed even when the model is
  # mathematically equivalent - preserving order keeps a seeded before/after
  # diff exactly reproducible rather than merely statistically similar.
  paste0("model {\n", paste(rev(blocks), collapse = ""), "  }")
}

BaHZING_Model <- function(formatted_data,
                          x,
                          covar = NULL,
                          exposure_standardization = NULL,
                          n.chains = 3,
                          n.adapt = 100,
                          n.iter.burnin = 1000,
                          n.iter.sample = 5000,
                          counterfactual_profiles = NULL,
                          q = NULL,
                          verbose = TRUE,
                          return_all_estimates = FALSE,
                          seed = NULL,
                          parallel = FALSE,
                          n.cores = NULL,
                          ROPE_range = c(-0.1, 0.1),
                          return_posterior = FALSE,
                          offset = TRUE) {

  if (!is.numeric(n.chains) || length(n.chains) != 1L || is.na(n.chains) ||
      !is.finite(n.chains) || n.chains < 1 || n.chains != floor(n.chains)) {
    stop("n.chains must be a positive integer.")
  }
  n.chains <- as.integer(n.chains)
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
       !is.finite(seed) || seed < 1 || seed != floor(seed) ||
       seed > .Machine$integer.max - n.chains)) {
    stop("seed must be a positive integer small enough to create one seed per chain.")
  }
  if (!is.null(n.cores) &&
      (!is.numeric(n.cores) || length(n.cores) != 1L || is.na(n.cores) ||
       !is.finite(n.cores) || n.cores < 1 || n.cores != floor(n.cores))) {
    stop("n.cores must be NULL or a positive integer.")
  }
  if (!is.logical(offset) || length(offset) != 1L || is.na(offset)) {
    stop("offset must be either TRUE or FALSE.")
  }

  # 1. Check input data ----
  # Extract metadata file from formatted data
  exposure_covar_dat <- data.frame(formatted_data$Table)

  # Create covariate dataframe
  if (!is.null(covar)){
    W <- data.frame(exposure_covar_dat[covar])
    Q <- ncol(W)
  }

  # Create exposure dataframe
  if(!all(x %in% colnames(exposure_covar_dat))) {
    stop("Not all exposured are found in the formatted data")
  }
  X <- exposure_covar_dat[x]
  P <- ncol(X)

  # Set exposure_standardization if missing
  if(is.null(exposure_standardization)) {
    exposure_standardization = "standard_normal"
  }
  # Check exposure_standardization is within the bounds
  if(!(exposure_standardization %in% c("standard_normal", "quantile", "none"))){
    stop("exposure_standardization must be either standard_normal, quantile, or none")
  }

  # Give warning if exposure_standardization == "standard_normal" and q is specified
  if(exposure_standardization == "standard_normal" & !is.null(q)){
    message("Note: q is not required when exposure_standardization is standard_normal. q will be ignored.")
  }

  # If standardization is "none", then check to make sure that counterfactual_profiles is specified
  if(exposure_standardization  == "none" & is.null(counterfactual_profiles)){
    stop("counterfactual_profiles must be speficied if exposure_standardization is none.")
  }

  # 2. Counterfactual profiles -----------------------------------
  ## 2.a Set counterfactual_profiles if missing --------
  if(is.null(counterfactual_profiles)) {
    if(exposure_standardization == "standard_normal") {
      counterfactual_profiles <- c(-0.5,0.5)
    } else {
      counterfactual_profiles <- c(0,1)
    }
  }

  ## 2.b Check counterfactual profiles structure ----
  if (is.matrix(counterfactual_profiles)) { # Checks if matrix
    if (nrow(counterfactual_profiles) != 2) {
      stop("counterfactual_profiles must have 2 rows when provided as a matrix.")
    }
    if (ncol(counterfactual_profiles) != P) {
      stop(paste0("When provided as a matrix, the number of columns in counterfactual_profiles must be equal to the number of exposures in the model."))
    }
    if (!is.numeric(counterfactual_profiles)) {
      stop("counterfactual_profiles must be numeric.")
    }
  } else { # Checks if is vector
    if (is.vector(counterfactual_profiles)) {
      if (!is.numeric(counterfactual_profiles)) {
        stop("counterfactual_profiles must be numeric.")
      }
      if (length(counterfactual_profiles) != 2) {
        stop(paste0("counterfactual_profiles must have 2 elements when provided as a vector."))
      }
    } else {
      stop("counterfactual_profiles must be either a numeric 2xP matrix or a numeric vector with length P.")
    }
  }

  ## 2.c Set profiles --------------------------------
  if(is.matrix(counterfactual_profiles)) {
    profiles = counterfactual_profiles
  } else {
    profiles <- rbind(rep(counterfactual_profiles[1], P),
                      rep(counterfactual_profiles[2], P))
  }

  chain_inits <- if (!is.null(seed)) {
    lapply(seq_len(n.chains), function(i) {
      list(.RNG.name = "base::Wichmann-Hill",
           .RNG.seed = seed + i)
    })
  } else {
    NULL
  }

  if (isTRUE(parallel) && is.null(n.cores)) {
    detected_cores <- parallel::detectCores(logical = FALSE)
    if (is.na(detected_cores) || detected_cores < 1) {
      detected_cores <- 1
    }
    n.cores <- max(1, min(n.chains, detected_cores))
  }
  if (!is.null(n.cores)) {
    n.cores <- min(n.chains, as.integer(n.cores))
  }

  # Give warning if using quantiles but counterfactuals < 0
  if(exposure_standardization=="quantile" &
     any(counterfactual_profiles<0 |
         counterfactual_profiles > q+1)){
    warning("Note: Quantiles are used, but counterfactual_profiles includes values outside of range. Estimates will be calculated based on the specified counterfactual_profiles, but results may not be interpretable.")
  }


  # 3. Scale exposures --------------------------------
  # Set default q if not provided
  if(is.null(q) & exposure_standardization == "quantile") {
    q = 4
  }

  # If using quantiles, quantize X
  if(exposure_standardization=="quantile") {
    probs <- seq(0, 1, length.out = q + 1)
    X.q <- apply(X, 2, function(v) {
      cut(v, breaks = c(-Inf, quantile(v, probs = probs, include.lowest = FALSE)), labels = FALSE)
    })
  }

  #If not quantized and not standardized, scale X
  if(exposure_standardization=="standard_normal") {
    X.q <- scale(X)
  }

  #If none, no scaling
  if(exposure_standardization == "none") {
    X.q <- X
  }

  # 4. Format microbiome matricies ----
  taxa_levels <- formatted_data$taxa_levels
  if (is.null(taxa_levels)) {
    stop("formatted_data has no $taxa_levels; regenerate it with the current Format_BaHZING().")
  }
  validate_taxa_levels(taxa_levels)
  narrowest_level <- taxa_level_name(taxa_levels, length(taxa_levels))

  # Which columns of exposure_covar_dat are taxon-count (outcome) data,
  # read directly from Format_BaHZING()'s own authoritative record rather
  # than inferred via string-matching (e.g. assuming every such column
  # contains "k__", which silently breaks if Kingdom/Domain isn't part of
  # the input taxonomy at all - Format_BaHZING() already computes this
  # exact list when it builds these columns, so there's no need to guess).
  taxon_columns <- formatted_data$taxon_columns
  if (is.null(taxon_columns) || length(taxon_columns) == 0) {
    stop("formatted_data has no $taxon_columns; regenerate it with the current Format_BaHZING().")
  }

  #Create outcome dataframe
  Y <- exposure_covar_dat[, taxon_columns, drop = FALSE]
  N <- nrow(Y)
  R <- ncol(Y)

  if (offset) {
    library_size <- rowSums(Y)
    if (any(!is.finite(library_size) | library_size <= 0)) {
      stop("Library sizes must be finite and greater than zero when offset = TRUE.")
    }
    L <- matrix(library_size, ncol = 1L)
  }

  # One data.frame + taxon-count per level except the narrowest (species-role)
  # level, e.g. GenusData/Genus.R, FamilyData/Family.R, ... - read from the
  # binary incidence matrices Format_BaHZING() built, keyed generically via
  # hierarchy_matrix_name() so this never has to hardcode a level name.
  level_data <- list()
  for (i in seq_len(length(taxa_levels) - 1)) {
    lvl <- taxa_level_name(taxa_levels, i)
    mat_name <- hierarchy_matrix_name(taxa_levels, i)
    if (is.null(formatted_data[[mat_name]])) {
      stop(sprintf("formatted_data is missing '%s'; regenerate it with Format_BaHZING().", mat_name))
    }
    df <- as.data.frame(t(formatted_data[[mat_name]]))
    level_data[[lvl]] <- list(data = df, R = ncol(df))
  }

  # 5. Return "Sanity" Messages ----
  if(verbose == TRUE){
    message("#### Checking input data ####")
    message("Exposure and Covariate Data:")
    message(paste0("- Total sample size: ", N))
    message(paste0("- Number of exposures: ", P))

    message("Microbiome Data:")
    # Narrowest (of the non-species-role levels) to broadest, e.g.
    # genus, family, order, class, phylum for the default hierarchy -
    # rev() of level_data's broadest-first insertion order above.
    for (lvl in rev(names(level_data))) {
      message(paste0("- Number of unique ", tolower(lvl), " in data: ", level_data[[lvl]]$R))
    }

    message("#### Running BaHZING with the following parameters #### ")
    if (exposure_standardization == "standard_normal"){
      message("Exposure standardization: Standard Normal")
    }
    if (exposure_standardization == "none"){
      message("Exposure standardization: None")
    }
    if (exposure_standardization == "quantile"){
      message(paste0("Exposure standardization: Quantiles, with q = ", q))
    }
    message("Library size offset: ", if (offset) "Included" else "Not included")
  }

  # 6. Run Model ----
  BHRM.microbiome <- .bahzing_build_model_text(
    taxa_levels,
    has_covar = !is.null(covar),
    has_offset = offset
  )

  ### Run JAGs Estimation ----
  # set up for JAGs based on taxonomy
  jdata <- list(N=N, Y=Y, R=R, X.q=X.q, P=P, profiles=profiles)
  if (offset) {
    jdata$L <- L
  }
  for (lvl in names(level_data)) {
    jdata[[paste0(lvl, "Data")]] <- level_data[[lvl]]$data
    jdata[[paste0(lvl, ".R")]] <- level_data[[lvl]]$R
  }
  if (!is.null(covar)) {
    jdata$Q <- Q
    jdata$W <- W
  }

  # Narrowest -> broadest prefix order, matching the original hardcoded var.s
  # (e.g. species, genus, family, order, class, phylum for the default
  # hierarchy): narrowest_level first, then level_data's names (broadest ->
  # narrowest insertion order, built in section 4) reversed.
  level_prefixes <- tolower(c(narrowest_level, rev(names(level_data))))
  var.s <- c(paste0(level_prefixes, ".beta"), paste0(level_prefixes, ".beta.zero"),
             paste0(level_prefixes, ".psi"), paste0(level_prefixes, ".psi.zero"),
             "disp")

  use_parallel <- isTRUE(parallel) && .Platform$OS.type != "windows" &&
    n.cores > 1L && n.chains > 1L

  if (use_parallel && is.null(chain_inits)) {
    # Forked processes inherit the parent's state. Give JAGS a distinct,
    # parent-generated RNG seed for each chain even when reproducibility was
    # not requested, so separate chains cannot accidentally start alike.
    chain_seeds <- sample.int(.Machine$integer.max - 1L, n.chains)
    chain_inits <- lapply(chain_seeds, function(s) {
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = s)
    })
  }

  if (use_parallel) {
    chain_list <- parallel::mclapply(seq_len(n.chains), function(chain_index) {
      chain_model <- jags.model(file=textConnection(BHRM.microbiome),
                                data=jdata,
                                n.chains=1,
                                n.adapt=n.adapt,
                                quiet=TRUE,
                                inits=chain_inits[[chain_index]])
      update(chain_model, n.iter=n.iter.burnin, progress.bar="none")
      chain <- coda.samples(model=chain_model,
                            variable.names=var.s,
                            n.iter=n.iter.sample,
                            thin=1,
                            progress.bar="none")
      chain[[1L]]
    }, mc.cores = n.cores)

    # mclapply records a failed fork as a try-error. Report the underlying
    # chain error directly instead of allowing a later, misleading
    # "Arguments must be mcmc objects" failure from coda.
    failed <- vapply(chain_list, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      messages <- vapply(chain_list[failed], function(x) {
        condition <- attr(x, "condition")
        if (is.null(condition)) as.character(x) else conditionMessage(condition)
      }, character(1))
      stop("Parallel chain(s) failed: ",
           paste(unique(messages), collapse = "; "))
    }
    model.fit <- coda::as.mcmc.list(chain_list)
  } else {
    model.fit <- jags.model(file=textConnection(BHRM.microbiome),
                            data=jdata,
                            n.chains=n.chains,
                            n.adapt=n.adapt,
                            quiet=F,
                            inits=chain_inits)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit,
                              variable.names=var.s,
                              n.iter=n.iter.sample,
                              thin=1,
                              progress.bar="text")
  }

  # 7. summarize results -------------------------------------------------------
  ## Calculate Mean, SD, and quantiles ----
  r <- summary(model.fit)
  results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

  ## Calculate HPD intervals (not included in current version) ------
  # x1 <- HPDinterval(model.fit[[1]], prob = 0.95)  %>% as.data.frame()

  ## Calculate p-values  ------
  # Operate on the native matrix (avoids coercing a large posterior matrix to
  # a data frame), and compute p_direction/p_rope/p_map in a single pass per
  # column instead of three separate apply() scans over the same columns.
  post_dist <- as.matrix(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]])), drop = FALSE]

  p_stats <- apply(post_dist, 2, function(x) {
    c(pdir  = as.numeric(p_direction(x, threshold = 0.05)),
      prope = p_rope(x, rope = ROPE_range)$p_ROPE,
      pmap  = as.numeric(p_map(x)))
  })

  # Get dataframe of p-values
  p_value_df <- data.frame(name = colnames(post_dist), t(p_stats), row.names = NULL)

  # # "Significant" example
  # bayestestR::describe_posterior(post_dist$`family.beta[1,1]`)
  # (pdirection <- bayestestR::p_direction(post_dist$`family.beta[1,1]`, as_p = TRUE)) %>% as.numeric()
  # (prope <- bayestestR::p_rope(post_dist$`family.beta[1,1]`, rope = c(-0.1, 0.1), as_p = TRUE) %>% as.numeric())
  # (pmap <- bayestestR::p_map(post_dist$`family.beta[1,1]`)) %>% as.numeric()
  # (ps <- bayestestR::p_significance(post_dist$`family.beta[1,1]`, rope = c(-0.1, 0.1), as_p = TRUE))
  #
  # # "Non-significant" example
  # bayestestR::describe_posterior(post_dist$`species.beta[189,2]`)
  # bayestestR::p_direction(post_dist$`species.beta[189,2]`, as_p = TRUE)
  # bayestestR::p_rope(post_dist$`species.beta[189,2]`, rope = c(-0.1, 0.1))

  ## Create "component" variable ----
  rn_results <- rownames(results)
  results <- results %>%
    mutate(
      component=case_when(
        grepl("zero",rn_results) ~ "Zero-inflation model coefficients",
        grepl("beta",rn_results) ~ "Count model coefficients",
        grepl("psi",rn_results)  ~ "Count model coefficients",
        grepl("disp",rn_results) ~ "Dispersion",
        #grepl("omega",rn_results) ~ "Omega",
        TRUE ~ "Other"))

  # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
  results <- results %>%
    # mutate(bci_inc_zero = case_when(
    #   grepl("disp",rownames(.)) ~ NA_character_,
    #   grepl("omega",rownames(.)) ~ NA_character_,
    #   (X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0) ~ "BCI ",
    #   TRUE ~ "N.S.")) %>%
    rename(estimate = Mean,
           bci_lcl = X2.5.,
           bci_ucl = X97.5.)

  # Combine results df with p-value df
  results$name = rownames(results)
  results2 <- left_join(results, p_value_df, by = "name")
  rownames(results2) <- results2$name
  results2 <- results2 %>% select(-name)
  # Rownames are fixed from here through the filter() below (no reordering
  # or row-adding operations occur in between), so cache once instead of
  # recomputing rownames(results2) on every grepl() call below.
  rn_results2 <- rownames(results2)

  # Calculate odds ratios-- removed --
  # results2 <- results2 %>%
  #   mutate(OR=exp(Mean),
  #          OR.ll=exp(X2.5.),
  #          OR.ul=exp(X97.5.))

  #Format output names
  # One taxon-name vector per level, keyed by level name (proper case, e.g.
  # "Genus") - narrowest level uses Y's colnames (the actual outcome/species
  # columns), every other level uses its level_data data.frame's colnames
  # (built in section 4).
  level_labels <- lapply(taxa_levels, function(lvl) {
    if (identical(lvl, narrowest_level)) colnames(Y) else colnames(level_data[[lvl]]$data)
  })
  names(level_labels) <- taxa_levels
  exposure <- colnames(X)
  results2$taxa_index <- str_remove(rn_results2,"..*\\[")
  results2$taxa_index <- str_remove(results2$taxa_index,",.*$")
  results2$taxa_index <- str_remove(results2$taxa_index,"]")
  results2$taxa_index <- as.numeric(results2$taxa_index)
  results2$Exposure.Index <- str_remove(rn_results2,"..*\\,")
  results2 <- results2 %>%
    mutate(Exposure.Index=ifelse(
      grepl("disp",Exposure.Index)|
        grepl("omega",Exposure.Index)|
        grepl("psi",Exposure.Index),
      NA,Exposure.Index))
  results2$Exposure.Index <- str_remove(results2$Exposure.Index,"]")
  results2$Exposure.Index <- as.numeric(results2$Exposure.Index)

  # Resolve which taxonomic level a JAGS parameter row belongs to, from its
  # rowname's "<level>."-prefixed parameter name (e.g. "genus.beta[12,1]" ->
  # "Genus") - anchored to the start of the parameter name so a level whose
  # name is a substring of another's can't silently mis-attribute (`stop()`s
  # instead). "disp"/"omega" rows have no level prefix at all - they're
  # dispersion/legacy parameters computed per narrowest-level taxon, so they
  # resolve to the narrowest level, matching the original model's behavior.
  resolve_level_for_row <- function(rowname) {
    if (grepl("disp", rowname) || grepl("omega", rowname)) {
      return(narrowest_level)
    }
    lower_levels <- tolower(taxa_levels)
    hits <- taxa_levels[vapply(lower_levels, function(p) {
      grepl(paste0("^", p, "\\."), rowname)
    }, logical(1))]
    if (length(hits) == 0) {
      stop(sprintf("Unrecognized taxonomic parameter '%s' (taxa_levels: %s).",
                   rowname, paste(taxa_levels, collapse = ", ")))
    }
    if (length(hits) > 1) {
      stop(sprintf("Ambiguous taxonomic parameter '%s' matches multiple levels: %s.",
                   rowname, paste(hits, collapse = ", ")))
    }
    hits
  }

  row_levels <- vapply(rownames(results2), resolve_level_for_row, character(1), USE.NAMES = FALSE)
  results2$domain <- row_levels
  results2$taxa_full <- mapply(function(lvl, idx) level_labels[[lvl]][idx],
                               row_levels, results2$taxa_index)
  results2$exposure <- paste0(exposure[results2$Exposure.Index])

  # Modify Exposure variable
  results2 <- results2 %>%
    mutate(exposure=case_when(
      grepl("psi",rn_results2) ~ "Mixture",
      grepl("disp",rn_results2) ~ "Dispersion",
      #grepl("omega",rn_results2) ~ "Omega",
      TRUE ~ exposure))

  # Get Taxa name from full taxonomic string
  results2 <- results2 %>%
    mutate(taxa_name = sub(".*__", "", taxa_full))

  # Remove "disp" and "omega" estimates
  if(!return_all_estimates){
    results2 <- results2 %>%
      filter(!grepl("disp",rn_results2),
             !grepl("omega",rn_results2))
  }

  # Remove rownames
  rownames(results2) <- NULL

  # Select final variables
  results2 <- results2 %>%
    select(taxa_full, taxa_name, domain, exposure,component,
           estimate,bci_lcl,bci_ucl,pdir,prope,pmap)

  if (return_posterior) {
    return(list(results = results2, posterior = model.fit))
  }
  return(results2)
}
