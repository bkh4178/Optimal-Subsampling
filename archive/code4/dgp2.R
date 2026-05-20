# dgp.R
# Data Generating Process — Simulation Section v5
#
# Scenarios:
#   S1 - Homoscedastic baseline        : sigma2 = 1
#   S2 - Variance-leverage aligned     : sigma2 = 1 + gamma * (ell / mean(ell))
#   S3 - Group-wise heteroskedasticity : sigma2 = 1 + 2 * G_s3  (G ~ Bernoulli(0.5), uncentered)
#
# 변경사항 (v4 → v5):
#   S2 DGP 수정 — 기존 exp(gamma * ell_score)에서 ell_score가 rank 기반 정규화라
#   gamma를 낮춰도 sigma2와 ell의 rank correlation이 항상 거의 1로 고정.
#   gamma가 실질적인 contrast 강도를 조절하도록 linear form으로 변경:
#     sigma2 = 1 + gamma * (ell / mean(ell))
#     - gamma = 0: homoskedastic
#     - gamma 커질수록 high-leverage 관측치의 분산이 linear하게 증가
#     - rank alignment 강하게 유지, exp 폭발 없음, 해석 단순
#   기존 ell_rank / ell_score 계산 제거.
#
# Settings: N = 10000, p = 9 (Normal 3 / Discrete Uniform 3 / Bernoulli 3)
# beta0 = rep(1, 9), gamma controls S2 leverage-variance contrast
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

  # Bernoulli: uncentered raw 따로 저장 (S3 variance DGP용)
  X_bern_raw <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
  X_bern     <- sweep(X_bern_raw, 2, colMeans(X_bern_raw), FUN = "-")

  X    <- cbind(X_norm, X_dunif, X_bern)
  G_s3 <- X_bern_raw[, 1]   # uncentered binary indicator: 0 또는 1
  # 열 순서: 1~3 Normal / 4~6 Discrete Uniform / 7~9 Bernoulli(centered)

  # Step 2. A_hat 및 leverage score
  A_hat     <- crossprod(X) / N
  A_hat_inv <- solve(A_hat)

  ell <- rowSums((X %*% A_hat_inv) * X)
  h   <- ell / N

  # Step 3. Scenario별 오차 분산
  # S2: gamma controls leverage-variance contrast (linear form)
  #   gamma=0 → homoskedastic, gamma 클수록 high-leverage 관측치 분산 증가
  # S3: G_s3=0 → sigma2_raw=1, G_s3=1 → sigma2_raw=3
  #   정규화 후: G=0 그룹 ≈ 0.5, G=1 그룹 ≈ 1.5 (약 3배 group-wise variance contrast)
  sigma2_raw <- switch(case,
    "1" = rep(1, N),
    "2" = 1 + gamma * (ell / mean(ell)),
    "3" = 1 + 2 * G_s3
  )

  sigma2 <- as.vector(sigma2_raw / mean(sigma2_raw))

  # Step 4. 반응변수 생성
  mu      <- as.vector(X %*% beta0)
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2))
  y       <- mu + epsilon

  list(
    X         = X,
    y         = y,
    mu        = mu,
    sigma2    = sigma2,
    ell       = ell,
    h         = h,
    A_hat     = A_hat,
    A_hat_inv = A_hat_inv,   # compute_Ainv_error diagnostic용
    beta0     = beta0,
    case      = case,
    gamma     = gamma,
    G_s3      = G_s3         # diagnostic용 (group별 sampling rate 확인 등)
  )
}