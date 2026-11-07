
PAPER_DIR   <- "/Users/annegiacobello/Desktop/EZ Paper "     
PROJECT_DIR <- file.path(PAPER_DIR, "sv hier 2 cond")
SCRIPT_04   <- file.path(PAPER_DIR, "04_compare_models_population.R")   

source(SCRIPT_04)   




PERSON <- list(

  EZ = list(
    map = c(nu_intercept = "nu_intercept", nu_c = "nu_c",
            alpha = "alpha_intercept", t0 = "tau_intercept"),
    transform = c(alpha = "exp")
  ),

  EZ_bias = list(
    map = c(nu_intercept = "nu_intercept", nu_c = "nu_c",
            alpha = "alpha_intercept", t0 = "tau_intercept", w = "w_person"),
    transform = c(alpha = "exp")
  ),

 
  DDM5 = list(
    map = c(nu_intercept = "nu_intercept", nu_c = "nu_c",
            alpha = "alpha_intercept", t0 = "t0", w = "w", sv = "sv"),
    transform = c(alpha = "exp")
  ),

  

  DDM6 = list(
    reader = read_summary("subj_summary.rds"),
    map = c(nu_intercept = "nu_intercept", nu_c = "nu_c", alpha = "a",
            t0 = "t0", w = "w", sv = "sv"),
    transform = c()
  ),

  DDM7 = list(
    reader = read_summary("subj_summary.rds"),
    map = c(nu_intercept = "nu_intercept", nu_c = "nu_c", alpha = "a",
            t0 = "t0", w = "w", sv = "sv"),
    transform = c()
  )
)

person_reader <- function(key) {
  if (!is.null(PERSON[[key]]$reader)) PERSON[[key]]$reader else MODELS[[key]]$reader
}


show_person_variables <- function(key, r = 1) {
  s <- person_reader(key)(MODELS[[key]]$dir, r)
  if (is.null(s)) { cat("no summary found for", key, "rep", r, "\n"); return(invisible(NULL)) }
  v <- s$variable[grepl("\\[\\d+\\]$", s$variable)]
  print(table(sub("\\[.*", "", v)))   # base name and number of persons
}




read_person_model <- function(key) {
  m <- MODELS[[key]]; pm <- PERSON[[key]]
  if (is.null(pm)) return(NULL)
  if (!dir.exists(m$dir)) { message(m$label, ": folder not found --> skipped"); return(NULL) }
  rdr <- person_reader(key)

  out <- bind_rows(lapply(setdiff(reps, m$exclude), function(r) {
    s <- rdr(m$dir, r)
    if (is.null(s)) return(NULL)
    s <- s[grepl("\\[\\d+\\]$", s$variable), ]
    s$base <- sub("\\[.*", "", s$variable)
    s$i    <- as.integer(sub(".*\\[(\\d+)\\]$", "\\1", s$variable))
    s <- s[s$base %in% pm$map, ]
   
    ids_file <- file.path(m$dir, sprintf("rep%02d", r), "subj_ids.rds")
    ids <- if (file.exists(ids_file)) readRDS(ids_file) else
      sort(readRDS(file.path(DATA_DIR, sprintf("rep%02d", r), "true_subject_parameters.rds"))$subj)
    s$subj  <- ids[s$i]
    s$param <- names(pm$map)[match(s$base, pm$map)]
    s$rep   <- r
    s
  }))
  if (nrow(out) == 0) { message(m$label, ": no person estimates found --> skipped"); return(NULL) }

  not_found <- setdiff(unname(pm$map), unique(out$base))
  if (length(not_found) > 0)
    warning(m$label, ": person vectors not found: ", paste(not_found, collapse = ", "),
            " --> check PERSON$", key, "$map (show_person_variables(\"", key, "\"))")

  for (p in names(pm$transform)) {
    f   <- match.fun(pm$transform[[p]])
    idx <- out$param == p
    for (cc in c("mean", "q5", "q95")) out[[cc]][idx] <- f(out[[cc]][idx])
  }
  out$model <- m$label
  out[, c("model", "rep", "subj", "param", "mean", "q5", "q95")]
}

est_p <- bind_rows(lapply(names(MODELS), read_person_model))
if (nrow(est_p) == 0) stop("no person estimates found")



truth_person_one <- function(r) {
  sp <- readRDS(file.path(DATA_DIR, sprintf("rep%02d", r), "true_subject_parameters.rds"))
  tibble(rep = r, subj = sp$subj,
         nu_intercept = sp$nu_intercept, nu_c = sp$nu_c,
         alpha = exp(sp$alpha_intercept), t0 = sp$t0,
         w = sp$w, sv = sp$sv, st0 = sp$st0) %>%
    pivot_longer(-c(rep, subj), names_to = "param", values_to = "truth")
}
truth_p <- bind_rows(lapply(reps, truth_person_one))




person_labels <- c(nu_intercept = "nu intercept", nu_c = "nu condition slope",
                   alpha = "alpha (boundary separation)", t0 = "t0 (seconds)",
                   w = "w (starting point)", sv = "sv (drift variability)",
                   st0 = "st0 (NDT variability, s)")

model_order_p <- intersect(sapply(MODELS, `[[`, "label"), unique(est_p$model))
colours_p     <- setNames(sapply(MODELS, `[[`, "colour"), sapply(MODELS, `[[`, "label"))[model_order_p]

pdat <- est_p %>%
  inner_join(truth_p, by = c("rep", "subj", "param")) %>%
  filter(param %in% names(person_labels)) %>%
  mutate(model = factor(model, levels = model_order_p),
         label = factor(person_labels[param], levels = person_labels)) %>%
  droplevels()


rec_rep <- pdat %>%
  group_by(label, model, rep) %>%
  summarise(r          = cor(truth, mean),
            slope      = coef(lm(mean ~ truth))[2],
            bias       = mean(mean - truth),
            rmse       = sqrt(mean((mean - truth)^2)),
            coverage90 = mean(truth >= q5 & truth <= q95),
            .groups = "drop")

r_lab <- rec_rep %>%
  group_by(label, model) %>%
  summarise(txt = sprintf("r = %.2f [%.2f, %.2f]", mean(r), min(r), max(r)), .groups = "drop") %>%
  group_by(label) %>%
  mutate(vj = -0.5 - 1.5 * (n() - rank(as.integer(model)))) %>%   # stack labels upward from bottom-right
  ungroup()


lims <- pdat %>%
  group_by(label) %>%
  summarise(lo = min(c(truth, mean)), hi = max(c(truth, mean)), .groups = "drop") %>%
  pivot_longer(c(lo, hi), values_to = "v")

legend_guide <- guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5, linetype = 0)))


axis_breaks <- list(scale_x_continuous(breaks = scales::breaks_extended(4)),
                    scale_y_continuous(breaks = scales::breaks_extended(4)))




p_ind_overlay <- ggplot(pdat, aes(x = truth, y = mean, colour = model)) +
  geom_blank(data = lims, aes(x = v, y = v), inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  geom_text(data = r_lab, aes(x = Inf, y = -Inf, label = txt, colour = model, vjust = vj),
            hjust = 1.05, size = 3, show.legend = FALSE, inherit.aes = FALSE) +
  facet_wrap(~ label, scales = "free", ncol = 4, labeller = label_wrap_gen(width = 30)) +
  axis_breaks +
  scale_colour_manual(values = colours_p, name = NULL) +
  legend_guide +
  labs(x = "True value", y = "Estimated value (posterior mean)",
       caption = sprintf("All persons pooled over reps. Lines: linear fit per model; red dashed: identity. r = mean [min, max] over reps.")) +
  theme_minimal() +
  theme(aspect.ratio = 1, legend.position = "bottom", strip.text = element_text(size = 10))


SHOW_IN_PANE <- FALSE
open_png <- function(f) if (Sys.info()[["sysname"]] == "Darwin") system2("open", shQuote(f))

if (SHOW_IN_PANE) print(p_ind_overlay)
n_rows_a <- ceiling(nlevels(pdat$label) / 4)
f_overlay <- file.path(PLOT_DIR, "individual_recovery_all_models_overlay.png")
ggsave(f_overlay, p_ind_overlay,
       width = 16, height = 4.2 * n_rows_a + 1.5, dpi = 300, bg = "white")
open_png(f_overlay)




combos   <- pdat %>% distinct(model, label)                      
lims_b   <- combos %>% left_join(lims, by = "label", relationship = "many-to-many")
r_lab_b  <- r_lab
n_params <- nlevels(pdat$label)
n_models <- nlevels(pdat$model)

p_ind_grid <- ggplot(pdat, aes(x = truth, y = mean, colour = model)) +
  geom_blank(data = lims_b, aes(x = v, y = v), inherit.aes = FALSE) +
  geom_abline(data = combos, aes(intercept = 0, slope = 1),
              colour = "red", linetype = "dashed") +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_text(data = r_lab_b, aes(x = Inf, y = -Inf, label = txt),
            hjust = 1.05, vjust = -0.6, size = 3, colour = "black", inherit.aes = FALSE) +
  facet_wrap(~ model + label, scales = "free", ncol = n_params, drop = FALSE,
             labeller = label_wrap_gen(width = 30, multi_line = TRUE)) +   # strip: model / parameter
  axis_breaks +
  scale_colour_manual(values = colours_p, name = NULL) +
  legend_guide +
  labs(x = "True value", y = "Estimated value (posterior mean)") +
  theme_minimal() +
  theme(aspect.ratio = 1, legend.position = "bottom", strip.text = element_text(size = 9))

if (SHOW_IN_PANE) print(p_ind_grid)
f_grid <- file.path(PLOT_DIR, "individual_recovery_all_models_grid.png")
ggsave(f_grid, p_ind_grid,
       width = 3.4 * n_params + 1, height = 3.6 * n_models + 1.2,
       dpi = 300, bg = "white", limitsize = FALSE)
open_png(f_grid)




rec_person_table <- rec_rep %>%
  group_by(label, model) %>%
  summarise(across(c(r, slope, bias, rmse, coverage90), mean), n_reps = n(), .groups = "drop") %>%
  arrange(label, model)

print(rec_person_table, n = Inf)
write.csv(rec_person_table, file.path(PLOT_DIR, "individual_recovery_all_models.csv"), row.names = FALSE)
