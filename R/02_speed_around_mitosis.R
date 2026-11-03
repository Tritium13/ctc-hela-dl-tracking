# =============================================================================
# 02_speed_around_mitosis.R
# Does nuclear speed change around cell division?
# Fluo-N2DL-HeLa, sequence 01; reference tracks vs. pipeline tracks (fine-tuned Cellpose-SAM + Trackastra).
#
# Input : results/steps_speed_seq01.csv   (notebook 07; one row per step)
# Output: results/r_speed_effects_seq01.csv
#         results/figures/r_speed_ratio_by_time.png
#         results/figures/r_speed_diagnostics.png
#
# Run from the repository root:
#   Rscript R/02_speed_around_mitosis.R | tee results/r_speed_output.txt
# Requires: lme4
# =============================================================================

suppressPackageStartupMessages(library(lme4))

WINDOW_FRAMES <- 16     # ±8 h around division (30 min per frame)
OFFSET        <- 0.5    # µm/h added before log, so that near-zero speeds do not dominate

s <- read.csv("results/steps_speed_seq01.csv")
s <- subset(s, !is.na(rel_frames) & abs(rel_frames) <= WINDOW_FRAMES & is.finite(speed_um_h))

# --- Time bins ----------------------------------------------------------------
# 1-h bins (2 frames) on both sides; the division step (rel_frames = 0) gets its own bin.
# Baseline = 4-8 h before division (bins -8 ... -5), pooled so that the reference level is not
# a single noisy hour. All other bins are compared with this baseline.
add_bins <- function(d) {
  d$time_bin <- ifelse(d$rel_frames < 0, -ceiling(-d$rel_frames / 2),
                ifelse(d$rel_frames > 0,  ceiling(d$rel_frames / 2), 0))
  d$bin <- ifelse(d$time_bin <= -5, "baseline", as.character(d$time_bin))
  d$bin <- factor(d$bin, levels = c("baseline", as.character(c(-4:-1, 0, 1:8))))
  d$log_speed <- log(d$speed_um_h + OFFSET)
  d
}

# --- Model ----------------------------------------------------------------------
# log(speed) ~ time bin, with two crossed random intercepts:
#   * division_id: all steps around one division share the cells' context (position, neighbours,
#     lineage) -> not independent;
#   * track: repeated steps of the same cell. A track can contribute to two divisions (born in one,
#     dividing in the next), so track is crossed with, not nested in, division.
# log scale: speeds are positive and right-skewed; effects become ratios (fold changes).
fit_source <- function(d) {
  m  <- lmer(log_speed ~ bin + (1 | division_id) + (1 | track), data = d, REML = TRUE)
  fe <- fixef(m)
  ci <- confint(m, parm = "beta_", method = "Wald")
  out <- data.frame(term = names(fe), est = unname(fe),
                    lo = ci[names(fe), 1], hi = ci[names(fe), 2], row.names = NULL)
  out <- out[out$term != "(Intercept)", ]
  out$bin_h    <- as.numeric(sub("^bin", "", out$term))
  out$ratio    <- exp(out$est)
  out$ratio_lo <- exp(out$lo)
  out$ratio_hi <- exp(out$hi)
  list(model = m, effects = out)
}

res <- list()
for (src in c("reference", "pipeline")) {
  d <- add_bins(subset(s, source == src))
  cat(sprintf("\n=== %s tracks ===\n", src))
  cat(sprintf("steps: %d | divisions: %d | tracks: %d | median speed: %.1f µm/h\n",
              nrow(d), length(unique(d$division_id)), length(unique(d$track)), median(d$speed_um_h)))
  f <- fit_source(d)
  vc <- as.data.frame(VarCorr(f$model))
  cat("variance components (log scale):\n")
  print(vc[, c("grp", "vcov")], row.names = FALSE, digits = 3)
  base <- exp(fixef(f$model)[["(Intercept)"]]) - OFFSET
  cat(sprintf("baseline speed (4-8 h before division, back-transformed): %.1f µm/h\n", base))
  cat("speed relative to baseline (ratio, 95% Wald CI):\n")
  tab <- f$effects[, c("bin_h", "ratio", "ratio_lo", "ratio_hi")]
  print(format(tab, digits = 3), row.names = FALSE)
  f$effects$source <- src
  f$data <- d
  res[[src]] <- f
}

# --- Robustness: does the pipeline reproduce the reference pattern? -----------
ref  <- res$reference$effects
pipe <- res$pipeline$effects
both <- merge(ref, pipe, by = "bin_h", suffixes = c("_ref", "_pipe"))
cat("\n=== Pipeline vs. reference ===\n")
cat(sprintf("correlation of bin effects (log scale): r = %.2f (n = %d bins)\n",
            cor(both$est_ref, both$est_pipe), nrow(both)))
cat(sprintf("bin with the highest speed: reference %+g h | pipeline %+g h\n",
            ref$bin_h[which.max(ref$est)], pipe$bin_h[which.max(pipe$est)]))
agree <- sign(both$lo_ref) == sign(both$lo_pipe) & sign(both$hi_ref) == sign(both$hi_pipe)
cat(sprintf("bins with the same conclusion (CI above / below / including 1): %d of %d\n",
            sum(agree), nrow(both)))

# --- Figure: ratio vs. time, both sources ------------------------------------------
png("results/figures/r_speed_ratio_by_time.png", width = 1500, height = 850, res = 150)
yl <- range(c(ref$ratio_lo, ref$ratio_hi, pipe$ratio_lo, pipe$ratio_hi, 1))
plot(NA, xlim = c(-4.5, 8.5), ylim = yl, log = "y",
     xlab = "time relative to division (h; 1-h bins, 0 = division step)",
     ylab = "speed relative to 4-8 h before division",
     main = "Nuclear speed around mitosis (sequence 01)")
abline(h = 1, lty = 2); abline(v = 0, lty = 3, col = "grey40")
for (k in seq_along(res)) {
  e  <- res[[k]]$effects
  dx <- c(-0.12, 0.12)[k]
  col <- c("#1f77b4", "#ff7f0e")[k]
  segments(e$bin_h + dx, e$ratio_lo, e$bin_h + dx, e$ratio_hi, col = col)
  points(e$bin_h + dx, e$ratio, pch = 19, col = col)
}
legend("topright", legend = names(res), col = c("#1f77b4", "#ff7f0e"), pch = 19, bty = "n")
invisible(dev.off())

# --- Diagnostics (reference model) -----------------------------------------------
m <- res$reference$model
png("results/figures/r_speed_diagnostics.png", width = 1600, height = 1200, res = 150)
par(mfrow = c(2, 2))
plot(fitted(m), resid(m), xlab = "fitted", ylab = "residual", main = "Residuals vs fitted (reference)")
abline(h = 0, lty = 2)
qqnorm(resid(m), main = "Normal Q-Q: residuals"); qqline(resid(m))
re <- ranef(m)$division_id[, 1]
qqnorm(re, main = "Normal Q-Q: division effects"); qqline(re)
hist(res$reference$data$log_speed, breaks = 40, main = "log(speed + 0.5)", xlab = "log(µm/h)")
invisible(dev.off())

# --- Save --------------------------------------------------------------------------
eff <- rbind(ref, pipe)[, c("source", "bin_h", "est", "lo", "hi", "ratio", "ratio_lo", "ratio_hi")]
write.csv(eff, "results/r_speed_effects_seq01.csv", row.names = FALSE)
cat("\nSaved: results/r_speed_effects_seq01.csv, results/figures/r_speed_*.png\n")
