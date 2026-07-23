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
# plugin 방법 추가
get_scores <- function(dat, sigma2_hat_param = NULL, sigma2_hat_nonparam = NULL) {
  N <- nrow(dat$X)
  scores <- list(
    "SRS"               = rep(1, N),
    "LEV-IPW"           = dat$ell,
    "OPT-homo"          = sqrt(dat$ell),
    "OPT-hetero-oracle" = sqrt(dat$sigma2_vec * dat$ell)
  )
  if (!is.null(sigma2_hat_param)) {
    scores[["OPT-hetero-plugin-param"]] <- sqrt(sigma2_hat_param * dat$ell)
  }
  if (!is.null(sigma2_hat_nonparam)) {
    scores[["OPT-hetero-plugin-nonparam"]] <- sqrt(sigma2_hat_nonparam * dat$ell)
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



# IBOSS 추가 (0714)
select_iboss <- function(X, k) {
  Xcov <- X[, -1, drop = FALSE]
  d    <- ncol(Xcov)
  n    <- nrow(Xcov)

  r <- suppressWarnings({
    if (k %/% d %/% 2 != k / d / 2) {
      message(sprintf("[IBOSS] k/d/2(k=%d,d=%d)가 정수가 아님 -> r=floor 사용", k, d))
    }
    max(k %/% d %/% 2, 1L)
  })

  get_idx_top_bottom <- function(z, r, exclude = integer(0)) {
    cand <- setdiff(seq_along(z), exclude)
    ord  <- cand[order(z[cand])]
    n_c  <- length(ord)
    r_j  <- min(r, n_c %/% 2L)
    c(ord[seq_len(r_j)], ord[(n_c - r_j + 1L):n_c])
  }

  idx_od <- get_idx_top_bottom(Xcov[, 1], r)
  for (j in 2:d) {
    idx_od <- c(idx_od, get_idx_top_bottom(Xcov[, j], r, exclude = idx_od))
  }
  idx_od
}

# IBOSS subdata unweighted OLS (deterministic, rep 무관)
sample_iboss <- function(dat, k) {
  idx <- select_iboss(dat$X, k)
  fit <- tryCatch(
    .lm.fit(x = dat$X[idx, , drop = FALSE], y = dat$y[idx]),
    error = function(e) NULL
  )
  beta_hat <- if (!is.null(fit)) {
    bhat <- as.numeric(fit$coefficients)
    if (length(bhat) == ncol(dat$X) && all(is.finite(bhat))) bhat else NULL
  } else NULL
  list(beta_hat = beta_hat, idx = idx, n_realized = length(idx))
}