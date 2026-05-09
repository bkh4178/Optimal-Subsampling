# dgp.R
# Data Generating Process — Simulation Section v3
#
# Scenarios:
#   S1 - Homoscedastic baseline          : sigma2 = 1
#   S2 - Variance-leverage aligned       : sigma2 = exp(gamma * ell_score)
#   S3 - Variance-leverage weakly sep.   : sigma2 = 1 + 2 * X[,7]^2
#
# Settings: N = 10000, p = 9 (Normal 3 / Discrete Uniform 3 / Bernoulli 3)
# beta0 = rep(1, 9), gamma = 0.7
# No intercept. All covariates mean-centered.

generate_data <- function(N = 10000, p = 9, case = 1, gamma = 0.7, beta0 = NULL) {

  if (p != 9) stop("This DGP is defined only for p = 9.")

  case <- as.character(case)
  if (!case %in% c("1", "2", "3")) stop("case must be 1, 2, or 3.")

  if (is.null(beta0)) beta0 <- rep(1, p)
  if (length(beta0) != p) stop("length(beta0) must equal p.")

  # Step 1. 공변량 생성 (3/3/3 구조)
  X_norm  <- matrix(rnorm(N * 3), N, 3)

  X_dunif <- matrix(sample(1:5, N * 3, replace = TRUE), N, 3)
  X_dunif <- sweep(X_dunif, 2, colMeans(X_dunif), FUN = "-")

  X_bern  <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
  X_bern  <- sweep(X_bern, 2, colMeans(X_bern), FUN = "-")

  X <- cbind(X_norm, X_dunif, X_bern)
  # 열 순서: 1~3 Normal / 4~6 Discrete Uniform / 7~9 Bernoulli

  # Step 2. A_hat 및 leverage score
  A_hat     <- crossprod(X) / N          # X'X / N = E[XX']의 sample analogue
  A_hat_inv <- solve(A_hat)

  ell <- rowSums((X %*% A_hat_inv) * X)  # ell_i = x_i' A^{-1} x_i
  h   <- ell / N

  # S2용 leverage score 안정화 (rank 기반 정규화, [-2, 2] clipping)
  ell_rank  <- rank(ell) / (N + 1)
  ell_score <- qnorm(ell_rank)
  ell_score <- pmax(pmin(ell_score, 2), -2)

  # Step 3. Scenario별 오차 분산
  sigma2_raw <- switch(case,
    "1" = rep(1, N),
    "2" = exp(gamma * ell_score),
    "3" = 1 + 2 * X[, 7]^2             # Bernoulli 열(7번) → leverage 기여 낮음
  )

  sigma2 <- as.vector(sigma2_raw / mean(sigma2_raw))  # 평균 noise level 통일

  # Step 4. 반응변수 생성
  mu      <- as.vector(X %*% beta0)
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2))
  y       <- mu + epsilon

  list(
    X      = X,
    y      = y,
    mu     = mu,
    sigma2 = sigma2,
    ell    = ell,
    h      = h,
    A_hat  = A_hat,
    beta0  = beta0,
    case   = case,
    gamma  = gamma
  )
}