# NOTE: the fixtures below are not a byte-for-byte copy of BaHZING's original
# hand-written model text - one deliberate change was made on top of it. The
# original always named the narrowest (species-role) level's precision nodes
# bare `tau`/`sigma` (every other level prefixes them, e.g. `genus.tau`) and
# its borrowed-mean node `mu.species` (named after itself, instead of
# `mu.genus`, the parent-naming convention every other level follows). Both
# are legacy inconsistencies with no functional weight - `tau`/`sigma`/`mu.*`
# are never monitored/extracted JAGS nodes - so they were made consistent
# with every other level's naming convention rather than preserved. The
# fixtures reflect that change; everything else here still matches the
# original model exactly.

# Normalizes JAGS text for structural comparison: strips "#..."-to-end-of-line
# comments (JAGS ignores these entirely - they carry no statistical meaning,
# and the original hand-written model text is inconsistent in its comment
# wording/capitalization from level to level, e.g. "#zero inflation" vs
# "#Zero inflation" vs "#zero infl" vs "#zero component") and strips all
# whitespace entirely, not just collapsing it - the original text is also
# inconsistent about e.g. spaces after commas (genus's block omits the space
# genus/family/order/class all otherwise use) - so two model strings that
# define the same JAGS model compare equal regardless of any cosmetic
# formatting differences.
normalize_jags <- function(text) {
  text <- gsub("#[^\n]*", "", text)
  text <- gsub("\\s+", "", text)
  text
}

read_fixture <- function(name) {
  paste(readLines(test_path("fixtures", name)), collapse = "\n")
}

test_that(".bahzing_build_model_text reproduces the no-offset with-covariates model", {
  generated <- .bahzing_build_model_text(
    default_taxa_levels, has_covar = TRUE, has_offset = FALSE
  )
  legacy <- read_fixture("legacy_model_with_covar.txt")
  expect_equal(normalize_jags(generated), normalize_jags(legacy))
})

test_that(".bahzing_build_model_text reproduces the no-offset model without covariates", {
  generated <- .bahzing_build_model_text(
    default_taxa_levels, has_covar = FALSE, has_offset = FALSE
  )
  legacy <- read_fixture("legacy_model_without_covar.txt")
  expect_equal(normalize_jags(generated), normalize_jags(legacy))
})

test_that(".bahzing_build_model_text adds the historical offset when requested", {
  generated <- .bahzing_build_model_text(
    default_taxa_levels, has_covar = TRUE, has_offset = TRUE
  )

  expect_match(
    generated,
    "log\\(lambda\\[i,r\\]\\).*log\\(L\\[i,1\\]\\)"
  )
  expect_match(
    generated,
    "logit\\(pi\\[i,r\\]\\).*log\\(L\\[i,1\\]\\)"
  )
})

test_that(".bahzing_build_model_text places the terminal (no-parent) block at the broadest level for a custom hierarchy", {
  custom_levels <- c("Phylum", "Family", "Genus", "ASV")
  generated <- .bahzing_build_model_text(custom_levels, has_covar = FALSE)

  # Broadest level (Phylum) should have no parent-derived mean - a flat
  # dnorm(0, ...) prior, and no mu node keyed by its own index (p.r).
  # ("mu.phylum[f.r,...]" legitimately appears in Family's block below,
  # since Family borrows from Phylum - only Phylum's own index (p.r) should
  # never appear as a mu subscript.)
  expect_true(grepl("phylum\\.beta\\[p\\.r,p\\] ~ dnorm\\(0, phylum\\.tau\\[p\\.r\\]\\)", generated))
  expect_false(grepl("mu\\.[a-z]+\\[p\\.r", generated))

  # Narrowest level (ASV, this hierarchy's data-likelihood role) should have
  # the actual data likelihood and borrow its prior mean from its parent,
  # Genus - not from a hardcoded "genus"/"Genus.R"/"GenusData". The borrowed-
  # mean node is named after the parent (mu.genus), same convention every
  # other level uses - not mu.asv (named after itself), which is what
  # BaHZING's original hand-written model did only for the narrowest level.
  expect_true(grepl("Y\\[i,r\\] ~ dnegbin", generated))
  expect_true(grepl("asv\\.beta\\[r,1:P\\]", generated))
  expect_true(grepl("mu\\.genus\\[r,p\\] <- inprod\\(genus\\.beta\\[1:Genus\\.R,p\\], GenusData\\[r,1:Genus\\.R\\]\\)", generated))

  # Mid-hierarchy level (Family) should borrow from its parent (Phylum)
  expect_true(grepl("mu\\.phylum\\[f\\.r,p\\] <- inprod\\(phylum\\.beta\\[1:Phylum\\.R,p\\], PhylumData\\[f\\.r,1:Phylum\\.R\\]\\)", generated))
})
