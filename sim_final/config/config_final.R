# config_final.R
# 전역 파라미터 정의 – sim_final 모든 스크립트에서 source() 로 읽어 사용
# 실행 위치: sim_final/ 폴더 기준

# ── 모집단 크기 ───────────────────────────────────────────────────────────────
N      <- 10000L
N_test <- 10000L

# ── Pilot 크기 (기본) ─────────────────────────────────────────────────────────
n0_default <- 200L

# ── 반복 횟수 ─────────────────────────────────────────────────────────────────
n_rep <- 500L

# ── Subsample 크기 ────────────────────────────────────────────────────────────
# Exp1: k 비교 실험
k_grid <- c(400L, 600L, 800L, 1000L, 2000L)

# Exp2: n0 sensitivity – k 고정
k_exp2  <- 1000L

n0_grid <- c(50L, 100L, 200L, 300L, 400L, 500L)

# ── DGP 파라미터 ──────────────────────────────────────────────────────────────
# X = cbind(1, Z1,...,Z9),  dim(X) = N x 10
# Z1~Z3 ~ N(0,1)
# Z4~Z6 ~ Bernoulli(0.5) - 0.5
# Z7~Z9 ~ Uniform(-1, 1)
# beta = (intercept, Z1, ..., Z9)
beta_true <- c(1, 1.5, 1.0, 0.75, 1.0, 0.7, 0.5, 0.5, 0.25, 0)

# 이분산 공식
dgp_sigma_label <- "sigma(Z) = 0.5 + |0.5*Z1 + 0.5*Z2 + Z6 + Z8|"

# ── Seed 전략 ─────────────────────────────────────────────────────────────────
# SEED_DATA : 전체 데이터 생성 (실험 내에서 고정)
# SEED_TEST : 테스트셋 생성 (실험 내에서 고정)
# SEED_REP  : rep r 의 seed = SEED_REP + r  (파일럿 + 서브샘플 결정)
SEED_DATA <- 123L
SEED_TEST <- 999L
SEED_REP  <- 1000L

# ── 저장 경로 (sim_final/ 기준) ───────────────────────────────────────────────
DIR_RAW  <- file.path("results", "raw")
DIR_EXP0 <- file.path(DIR_RAW, "exp0")
DIR_EXP1 <- file.path(DIR_RAW, "exp1")
DIR_EXP2 <- file.path(DIR_RAW, "exp2")
DIR_EXP3 <- file.path(DIR_RAW, "exp3")
DIR_EXP4 <- file.path(DIR_RAW, "exp4")

# ── 메서드 순서 ───────────────────────────────────────────────────────────────
METHOD_ORDER <- c("FULL","SRS", "LEV-IPW", "OPT-homo",
                  "OPT-hetero-oracle", "OPT-hetero-plugin")
