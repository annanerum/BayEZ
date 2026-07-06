# are the first 3-5 trials of each block slow?

library(dplyr)
library(ggplot2)
library(tidyr)

DATA_FILE  <- "krakow_data_standardized.csv"

PID_COL    <- "participant_index"
RT_COL     <- "rt"
COND_COL   <- "condition"          # congruent / incongruent
RESP_COL   <- "pre_acc"            # post-correct (1) / post-error (-1)

TRIAL_NUM_COL <- "trial_number"    # actual experiment-wide trial index (1-300)
BLOCK_SIZE    <- 60                # trials per block, block = ceiling(trial_number / BLOCK_SIZE)

N_FIRST_TRIALS <- 5                # how many "first trials" per block count as "early"
N_POSITIONS_TO_PLOT <- 15          # how many trial positions to show in the position plot

FILTER_IS_IN_SEQUENCE <- "is_in_sequence"
FILTER_RT_TOO_LONG    <- "rt_greater_than_1"
FILTER_LOG_RT_OUTLIER <- "log_rt_exceed_threshold"

OUT <- "output_rt_checks"
if (!dir.exists(OUT)) dir.create(OUT)


dark <- theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "#1e1e1e", colour = NA),
    panel.background  = element_rect(fill = "#1e1e1e", colour = NA),
    panel.grid.major  = element_line(colour = "#333333"),
    panel.grid.minor  = element_blank(),
    text              = element_text(colour = "white"),
    axis.text         = element_text(colour = "white"),
    strip.text        = element_text(colour = "white"),
    legend.background = element_rect(fill = "#1e1e1e"),
    legend.text       = element_text(colour = "white"),
    axis.title        = element_text(colour = "white")
  )

sv <- function(p, nm, w = 10, h = 7)
  ggsave(file.path(OUT, nm), p, width = w, height = h, dpi = 300, bg = "#1e1e1e")


cat("Loading data \n")
dat <- read.csv(DATA_FILE)

required <- c(PID_COL, RT_COL, COND_COL, TRIAL_NUM_COL, FILTER_IS_IN_SEQUENCE,
              FILTER_RT_TOO_LONG, FILTER_LOG_RT_OUTLIER)
missing  <- setdiff(required, names(dat))
if (length(missing) > 0) {
  stop("Missing column(s): ", paste(missing, collapse = ", "),
       "\nColumns found in csv: ", paste(names(dat), collapse = ", "))
}


dat <- dat %>%
  mutate(block     = ceiling(.data[[TRIAL_NUM_COL]] / BLOCK_SIZE),
         trial_pos = ((.data[[TRIAL_NUM_COL]] - 1) %% BLOCK_SIZE) + 1)

n_blocks_check <- dat %>% count(participant_index, block) %>% pull(n)
cat(sprintf("from %d to %d trials (expected up to %d)\n",
            min(n_blocks_check), max(n_blocks_check), BLOCK_SIZE))



dat$excluded <- !(dat[[FILTER_IS_IN_SEQUENCE]] == "True" &
                   dat[[FILTER_RT_TOO_LONG]]    == "False" &
                   dat[[FILTER_LOG_RT_OUTLIER]] == "False")

dat_clean <- dat[!dat$excluded, ]

cat(sprintf("trials total: %d , excluded: %d (%.1f%%) , remaining: %d\n",
            nrow(dat), sum(dat$excluded), 100 * mean(dat$excluded), nrow(dat_clean)))




p_hist_raw_vs_clean <- ggplot() +
  geom_histogram(data = dat, aes(x = .data[[RT_COL]], fill = "before exclusion"),
                  bins = 60, alpha = 0.5, position = "identity") +
  geom_histogram(data = dat_clean, aes(x = .data[[RT_COL]], fill = "after exclusion"),
                  bins = 60, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("before exclusion" = "steelblue",
                               "after exclusion"  = "goldenrod")) +
  labs(title = "RT distribution before vs. after exclusion filters",
       x = "RT (s)", y = "count", fill = NULL) +
  dark
sv(p_hist_raw_vs_clean, "hist_rt_before_after_exclusion.png")


if (COND_COL %in% names(dat_clean)) {
  p_hist_by_cond <- ggplot(dat_clean, aes(x = .data[[RT_COL]], fill = factor(.data[[COND_COL]]))) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    labs(title = "RT distribution by condition",
         x = "RT (s)", y = "count", fill = COND_COL) +
    dark
  sv(p_hist_by_cond, "hist_rt_by_condition.png")
}

if (all(c(COND_COL, RESP_COL) %in% names(dat_clean))) {
  p_hist_by_cell <- ggplot(dat_clean,
                            aes(x = .data[[RT_COL]],
                                fill = interaction(.data[[COND_COL]], .data[[RESP_COL]]))) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
    facet_wrap(~ interaction(.data[[COND_COL]], .data[[RESP_COL]])) +
    labs(title = "RT distribution by cell (condition x response-type)",
         x = "RT (s)", y = "count") +
    dark + theme(legend.position = "none")
  sv(p_hist_by_cell, "hist_rt_by_cell.png")
}


# FIRST TRIALS PER BLOCK

{
  dat_clean$is_early <- dat_clean$trial_pos <= N_FIRST_TRIALS

  # mean RT by trial position 
  pos_summary <- dat_clean %>%
    filter(trial_pos <= N_POSITIONS_TO_PLOT) %>%
    group_by(trial_pos) %>%
    summarise(mean_rt = mean(.data[[RT_COL]], na.rm = TRUE),
              se_rt   = sd(.data[[RT_COL]], na.rm = TRUE) / sqrt(n()),
              n       = n(), .groups = "drop")

  p_pos <- ggplot(pos_summary, aes(x = trial_pos, y = mean_rt)) +
    geom_ribbon(aes(ymin = mean_rt - se_rt, ymax = mean_rt + se_rt),
                fill = "steelblue", alpha = 0.3) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(colour = "goldenrod", size = 2) +
    geom_vline(xintercept = N_FIRST_TRIALS + 0.5, colour = "white", linetype = "dashed") +
    labs(title = "Mean RT by trial position within block (block = 60 trials)",
         subtitle = paste0("Dashed line = end of 'first ", N_FIRST_TRIALS, " trials' window"),
         x = "Trial position in block", y = "Mean RT (s), \u00b1 SE") +
    dark
  sv(p_pos, "rt_by_trial_position.png")

  # early versus rest boxplot (overall and by condition)
  p_box <- ggplot(dat_clean, aes(x = is_early, y = .data[[RT_COL]], fill = is_early)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.15) +
    scale_x_discrete(labels = c(`FALSE` = "rest of block",
                                 `TRUE`  = paste0("first ", N_FIRST_TRIALS, " trials"))) +
    scale_fill_manual(values = c(`TRUE` = "goldenrod", `FALSE` = "steelblue"), guide = "none") +
    labs(title = paste0("RT: first ", N_FIRST_TRIALS, " trials vs. rest of block"),
         x = NULL, y = "RT (s)") +
    dark
  sv(p_box, "rt_first_vs_rest.png")

  if (COND_COL %in% names(dat_clean)) {
    p_box_cond <- p_box + facet_wrap(vars(.data[[COND_COL]]))
    sv(p_box_cond, "rt_first_vs_rest_by_condition.png")
  }

  # descriptive numbers &  welch t-test
  desc <- dat_clean %>%
    group_by(is_early) %>%
    summarise(mean_rt = mean(.data[[RT_COL]], na.rm = TRUE),
              sd_rt   = sd(.data[[RT_COL]], na.rm = TRUE),
              n       = n(), .groups = "drop")
  cat("\nDescriptives (after exclusion):\n")
  print(desc)

  tt <- t.test(dat_clean[[RT_COL]] ~ dat_clean$is_early)
  cat(sprintf("\nWelch t-test, first %d trials vs. rest: t = %.2f, df = %.1f, p = %.4g\n",
              N_FIRST_TRIALS, tt$statistic, tt$parameter, tt$p.value))

  write.csv(desc, file.path(OUT, "first_trials_descriptives.csv"), row.names = FALSE)

}

cat("\n outputs saved to:", OUT, "\n")
cat("files: \n")
cat(paste0("  ", list.files(OUT), collapse = "\n"), "\n")
