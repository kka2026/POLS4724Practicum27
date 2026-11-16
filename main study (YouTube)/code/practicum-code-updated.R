# ---------------------------------------------
# MAIN STUDY
# ---------------------------------------------
library(randomizr)
library(estimatr)
library(ri2)
library(dplyr)
library(ggplot2)

rm(list = ls())
set.seed(111626)

topics <- c("morning routine", "night routine", "meal prep", "cleaning routine",
            "grocery haul", "sourdough bread", "baking from scratch",
            "skincare routine", "makeup tutorial", "hair care routine",
            "outfit ideas", "home decor", "garden tour", "homeschooling",
            "parenting tips", "toddler activities", "pregnancy tips",
            "marriage advice", "dating advice", "wedding planning",
            "budgeting tips", "laundry routine", "home organization",
            "easy dinner recipes", "sewing for beginners",
            "knitting for beginners", "workout routine", "self care routine",
            "day in the life", "kitchen tour")

n_units <- 30
n_treat <- 15
n_top   <- 20     # outcome counted among the top 20 results
n_sims  <- 10000  # randomization inference simulations
frame   <- "tradwife"

# random assignment (complete RA, 15/15) and random search order
Z <- complete_ra(N = n_units, m = n_treat)
assignments <- data.frame(
  topic_id = 1:n_units,
  topic    = topics,
  Z        = Z,
  query    = ifelse(Z == 1, paste(frame, topics), topics),
  order    = sample(1:n_units)
) |>
  arrange(order)
write.csv(assignments, "search_assignments.csv", row.names = FALSE)

# data: search_results.csv, one row per result (20 per topic):
#   topic_id, rank (1-20), video_id, title, channel, snippet,
#   political (0/1), note, searched_at
if (file.exists("search_results.csv")) {
  results <- read.csv("search_results.csv")
} else {
  # SAMPLE DATA so the code runs before data collection; delete afterwards
  results <- expand.grid(rank = 1:n_top, topic_id = 1:n_units) |>
    left_join(select(assignments, topic_id, Z), by = "topic_id") |>
    mutate(political = rbinom(n(), 1, 0.03 + 0.07 * Z)) |>
    select(topic_id, rank, political)
}

dat <- results |>
  filter(rank <= n_top) |>
  group_by(topic_id) |>
  summarise(Y = sum(political), n_results = n(), .groups = "drop") |>
  left_join(assignments, by = "topic_id") |>
  mutate(condition = factor(Z, levels = 0:1,
                            labels = c("Neutral query", "Conservative-coded query")))

# distribution of individual outcomes by condition
plot_dist <- ggplot(dat, aes(x = Y)) +
  geom_bar() +
  facet_wrap(~ condition, ncol = 1) +
  scale_x_continuous(breaks = 0:n_top, limits = c(-0.5, n_top + 0.5)) +
  theme_bw() +
  labs(x = "Explicitly political videos among top 20 results",
       y = "Number of topics")

# ATE, Neyman SE, 95% CI
fit <- difference_in_means(Y ~ Z, data = dat)

# randomization inference: sharp null of no effect, one-tailed (upper)
declaration <- declare_ra(N = n_units, m = n_treat)
ri_h1 <- conduct_ri(Y ~ Z, declaration = declaration, assignment = "Z",
                    sharp_hypothesis = 0, sims = n_sims, p = "upper",
                    data = dat)

results_table <- data.frame(
  estimate  = fit$coefficients, std.error = fit$std.error,
  conf.low  = fit$conf.low,     conf.high = fit$conf.high,
  ri_p_one_tailed = summary(ri_h1)$upper_p_value
)

# graphical depiction in the style of Coppock (2021)
fit_means <- lm_robust(Y ~ 0 + condition, data = dat)
means_df <- data.frame(
  condition = factor(levels(dat$condition), levels = levels(dat$condition)),
  Y = unname(fit_means$coefficients),
  conf.low = unname(fit_means$conf.low), conf.high = unname(fit_means$conf.high)
)
plot_coppock <- ggplot(means_df, aes(condition, Y)) +
  geom_point(data = dat, position = position_jitter(width = 0.15, height = 0.1),
             alpha = 0.4) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_y_continuous(limits = c(-0.5, n_top), breaks = seq(0, n_top, 2)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Explicitly political videos\n(of top 20 YouTube results)")

print(plot_dist); print(results_table); print(plot_coppock)


# ---------------------------------------------
# PILOT STUDY
# ---------------------------------------------
library(readxl)
library(tidyr)

p_assign <- read_excel("practicum_data.xlsx", sheet = "assignments") |>
  select(account, Z1, Z2) |> mutate(rounds = Z1 + Z2)
p_out    <- read_excel("practicum_data.xlsx", sheet = "outcomes")
p_videos <- read_excel("practicum_data.xlsx", sheet = "videos")

pilot <- p_assign |> left_join(p_out, by = "account")
pilot_complete <- pilot |>
  filter(!is.na(pol_wave0), !is.na(pol_wave1), !is.na(pol_wave2))

# each account's counts across waves (accounts 1-2)
pilot_traj <- pilot_complete |>
  mutate(change_0_1 = pol_wave1 - pol_wave0, change_1_2 = pol_wave2 - pol_wave1)

# baseline counts (accounts 1-3)
pilot_baseline <- pilot |> filter(!is.na(pol_wave0)) |> select(account, pol_wave0)

# liked-video descriptives by account, round and condition
pilot_videos <- p_videos |>
  filter(account %in% pilot_complete$account, !is.na(likes)) |>
  left_join(p_assign, by = "account") |>
  mutate(condition = if_else(if_else(round == 1, Z1, Z2) == 1,
                             "Conservative-coded", "Apolitical (cockapoo)")) |>
  group_by(account, round, condition) |>
  summarise(n_videos = n(),
            mean_likes = mean(likes), mean_comments = mean(comments),
            mean_followers = mean(followers),
            log_mean_likes = log(mean_likes + 1),
            log_mean_comments = log(mean_comments + 1),
            log_mean_followers = log(mean_followers + 1), .groups = "drop")

# trajectory plot
plot_pilot <- pilot_complete |>
  pivot_longer(starts_with("pol_wave"), names_to = "wave",
               names_prefix = "pol_wave", values_to = "political") |>
  mutate(account = factor(paste("Account", account, "| rounds treated =", rounds))) |>
  ggplot(aes(wave, political, group = account, linetype = account)) +
  geom_line() + geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 50)) + theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  labs(x = "Wave (0 = baseline, 1 = after round 1, 2 = after round 2)",
       y = "Explicitly political videos\n(of first 50 'For You' videos)",
       linetype = NULL)