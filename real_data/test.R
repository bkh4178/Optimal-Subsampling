install.packages('arrow')
library(arrow)

df <- read_parquet("yellow_tripdata_2023-01.parquet")

# 2023년 1월 yellow taxi 다운로드
url <- "https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2023-01.parquet"
download.file(url, destfile = "yellow_tripdata_2023-01.parquet", mode = "wb")

# 읽기
df <- read_parquet("yellow_tripdata_2023-01.parquet")
dim(df)
names(df)

# 기본 이상치 확인
summary(df$fare_amount)
summary(df$trip_distance)
table(df$RatecodeID)


summary(df_clean$fare_amount)
summary(df_clean$trip_distance)


df_clean <- df |>
  filter(RatecodeID == 1) |>
  filter(fare_amount > 2.5, fare_amount < 200) |>
  filter(trip_distance > 0, trip_distance < 50) |>
  collect() |>  # R 메모리로 가져오기
  mutate(
    trip_duration = as.numeric(difftime(
      tpep_dropoff_datetime,
      tpep_pickup_datetime,
      units = "mins")),
    hour = lubridate::hour(tpep_pickup_datetime)
  ) |>
  filter(trip_duration > 1, trip_duration < 120)

dim(df_clean)


# 1단계만
df2 <- df |>
  filter(RatecodeID == 1) |>
  collect()
dim(df2)



library(arrow)
library(dplyr)

df2 <- df |>
  filter(RatecodeID == 1) |>
  collect()
dim(df2)






library(arrow)
library(dplyr)
library(lubridate)

df_clean <- df |>
  filter(RatecodeID == 1) |>
  filter(fare_amount > 2.5, fare_amount < 200) |>
  filter(trip_distance > 0, trip_distance < 50) |>
  collect() |>
  mutate(
    trip_duration = as.numeric(difftime(
      tpep_dropoff_datetime,
      tpep_pickup_datetime,
      units = "mins")),
    hour = hour(tpep_pickup_datetime)
  ) |>
  filter(trip_duration > 1, trip_duration < 120)

dim(df_clean)


# 간단한 OLS 적합
fit <- lm(fare_amount ~ trip_distance + passenger_count + hour, data = df_clean)

# Residual plot (샘플링해서 그리기)
set.seed(42)
idx <- sample(nrow(df_clean), 10000)
plot(fitted(fit)[idx], residuals(fit)[idx],
     xlab = "Fitted", ylab = "Residuals",
     main = "Residuals vs Fitted", pch = 16, cex = 0.3, col = rgb(0,0,0,0.3))
abline(h = 0, col = "red")



library(lmtest)
install.packages('lmtest')
bptest(fit)
