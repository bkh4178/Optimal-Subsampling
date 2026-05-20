# 02_sampling_methods.R
# 각 method 의 sampling score 정의, Bernoulli 서브샘플 수행
#
# Score 정의:
#   SRS               : rep(1, N)         – 균등
#   LEV-IPW           : ell               – leverage 비례
#   OPT-homo          : sqrt(ell)         – 등분산 optimal
#   OPT-hetero-oracle : sqrt(sigma2 * ell)– 이분산 oracle
#   OPT-hetero-plugin : sqrt(sigma2_hat * ell) – plugin 추정치 사용

# 반환: method 명 → score 벡터의 named list
# sigma2_hat: plugin 추정치 (NULL 이면 OPT-hetero-plugin 제외)
get_scores <- function(dat, sigma2_hat = NULL) {
  N <- nrow(dat$X)
  scores <- list(
    "SRS"               = rep(1, N),
    "LEV-IPW"           = dat$ell,
    "OPT-homo"          = sqrt(dat$ell),
    "OPT-hetero-oracle" = sqrt(dat$sigma2_vec * dat$ell)
  )
  if (!is.null(sigma2_hat)) {
    scores[["OPT-hetero-plugin"]] <- sqrt(sigma2_hat * dat$ell)
  }
  scores
}

# score → pi → Bernoulli sampling
# 반환: list(idx, pi)
sample_bernoulli <- function(score, k) {
  pi  <- compute_pi(score, k)
  idx <- which(rbinom(length(pi), 1L, pi) == 1L)
  list(idx = idx, pi = pi)
}

# 전체 데이터 unweighted OLS (subsampling 없음)
# 반환: list(beta_hat, n_realized, n_clip)
sample_full <- function(dat) {
  fit <- tryCatch(
    .lm.fit(x = dat$X, y = dat$y),
    error = function(e) NULL
  )
  beta_hat <- if (!is.null(fit)) {
    bhat <- as.numeric(fit$coefficients)
    if (length(bhat) == ncol(dat$X) && all(is.finite(bhat))) bhat else NULL
  } else NULL
  list(beta_hat = beta_hat, n_realized = nrow(dat$X), n_clip = 0L)
}
