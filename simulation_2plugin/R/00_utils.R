# 00_utils.R
# 공통 유틸리티: compute_pi, fit_ipw, safe_cor, ensure_dir

# Bernoulli 서브샘플링 확률 계산
# score   : 각 관측치의 샘플링 score (≥ 0)
# k       : 목표 subsample 크기 (sum(pi) ≈ k)
# max_iter: pi > 1 클리핑 반복 상한
compute_pi <- function(score, k, max_iter = 100L) {
  N <- length(score)
  if (k <= 0) return(rep(0, N))
  if (all(!is.finite(score)) || sum(score, na.rm = TRUE) <= 0) {
    return(rep(min(k / N, 1), N))
  }
  score[!is.finite(score)] <- 0
  score <- pmax(score, 0)
  if (sum(score) <= 0) return(rep(min(k / N, 1), N))

  pi    <- score / sum(score) * k
  fixed <- rep(FALSE, N)

  for (iter in seq_len(max_iter)) {
    new_clip <- !fixed & (pi > 1)
    if (!any(new_clip)) break
    fixed[new_clip] <- TRUE
    pi[fixed] <- 1
    budget <- k - sum(fixed)
    free   <- !fixed
    if (budget <= 0 || !any(free)) { pi[free] <- 0; break }
    s_free <- sum(score[free])
    if (s_free <= 0) pi[free] <- budget / sum(free)
    else             pi[free] <- score[free] / s_free * budget
  }
  pmax(pmin(pi, 1), 0)
}

# IPW 추정량 (weight = 1/pi_i)
# 반환: 길이 ncol(X) 의 numeric 벡터, 실패 시 NULL
fit_ipw <- function(X, y, idx, pi) {
  if (length(idx) < ncol(X)) return(NULL)
  w <- 1 / pmax(pi[idx], 1e-10)
  if (any(!is.finite(w)) || any(w <= 0)) return(NULL)
  fit <- tryCatch(
    lm.wfit(x = X[idx, , drop = FALSE], y = y[idx], w = w),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  bhat <- as.numeric(fit$coefficients)
  if (length(bhat) != ncol(X) || any(!is.finite(bhat))) return(NULL)
  bhat
}

# 피어슨 상관: NA-안전, 분산 0 방어
safe_cor <- function(a, b, method = "pearson") {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  cor(a[ok], b[ok], method = method)
}

# 폴더 생성 (이미 있으면 무시)
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

  invisible(path)
}
