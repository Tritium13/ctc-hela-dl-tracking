# =============================================================================
# 01_paired_comparison.R
# Paired comparison of pretrained vs. fine-tuned Cellpose-SAM, per annotated nucleus
# (Fluo-N2DL-HeLa, sequence 01, gold truth).
#
# Input : results/per_nucleus_iou_seq01.csv   (notebook 05; one row per nucleus and model)
# Output: results/r_paired_comparison_summary.csv
#         results/figures/r_diagnostics.png
#         results/figures/r_iou_change_by_division.png
#
# Run from the repository root:
#   Rscript R/01_paired_comparison.R | tee results/r_paired_comparison_output.txt
#
# Requires: lme4 (everything else is base R)
# =============================================================================

suppressPackageStartupMessages(library(lme4))

# --- 1. Data -----------------------------------------------------------------
# Each annotated nucleus appears twice (once per model) -> paired design.
d <- read.csv("results/per_nucleus_iou_seq01.csv")
d$nucleus    <- paste(d$frame, d$gt_label, sep = "_")   # unique nucleus ID
d$matched    <- as.logical(d$matched)                    # pandas writes True/False
d$full_frame <- as.logical(d$full_frame)

# Wide format: one row per nucleus, values of both models side by side.
# Covariates are identical for both models, so they are taken once.
covars <- unique(d[, c("nucleus", "frame", "full_frame", "area_px", "frames_to_division")])
pre <- d[d$model == "pretrained", c("nucleus", "best_iou", "matched")]
ft  <- d[d$model == "fine-tuned", c("nucleus", "best_iou", "matched")]
names(pre)[2:3] <- c("iou_pre", "matched_pre")
names(ft)[2:3]  <- c("iou_ft",  "matched_ft")
w <- merge(merge(covars, pre, by = "nucleus"), ft, by = "nucleus")
stopifnot(nrow(w) == length(unique(d$nucleus)))          # every nucleus exactly once

# Near a division: within 2 frames (= 1 h) of the cell's birth or its own division.
# Nuclei without a division event in the movie (NA) count as "not near".
w$near_division <- !is.na(w$frames_to_division) & w$frames_to_division <= 2

cat("Nuclei:", nrow(w), "| frames:", length(unique(w$frame)),
    "| near division:", sum(w$near_division), "\n")
cat(sprintf("Mean best IoU  pretrained %.3f | fine-tuned %.3f\n", mean(w$iou_pre), mean(w$iou_ft)))
cat(sprintf("Detected (IoU > 0.5)  pretrained %d | fine-tuned %d\n\n",
            sum(w$matched_pre), sum(w$matched_ft)))

# --- 2. Detection: exact McNemar test ------------------------------------------
# Paired yes/no outcome. Only discordant nuclei carry information:
#   lost   = found by pretrained, missed by fine-tuned
#   gained = missed by pretrained, found by fine-tuned
# Under H0 (no difference) a discordant nucleus is equally likely to be lost or gained
# -> exact binomial test with p = 0.5 on the discordant pairs.
# Limitation: nuclei of the same frame are treated as independent here.
cat("=== Detection (2x2 table) ===\n")
print(table(pretrained = w$matched_pre, fine_tuned = w$matched_ft))
lost   <- sum(w$matched_pre & !w$matched_ft)
gained <- sum(!w$matched_pre & w$matched_ft)
cat(sprintf("\nlost: %d | gained: %d\n", lost, gained))
if (lost + gained > 0) {
  mc <- binom.test(gained, lost + gained, p = 0.5)
  cat(sprintf("exact McNemar: p = %.3f (proportion gained among discordant = %.2f, 95%% CI %.2f-%.2f)\n",
              mc$p.value, mc$estimate, mc$conf.int[1], mc$conf.int[2]))
}
cat("\nDiscordant nuclei by distance to division:\n")
print(table(near_division = w$near_division,
            change = ifelse(w$matched_pre & !w$matched_ft, "lost",
                     ifelse(!w$matched_pre & w$matched_ft, "gained", "unchanged"))))

# --- 3. Outline quality: linear mixed model on paired differences --------------
# Only nuclei detected by BOTH models (IoU of a missed nucleus measures something else).
# Outcome: difference in logit(IoU), fine-tuned minus pretrained.
#   * Differencing removes nucleus-to-nucleus variation (the pairing).
#   * logit: IoU is bounded at 1 and values cluster near 0.85-0.9; on the logit scale
#     differences are comparable across the range and residuals closer to normal.
#   * Random intercept for frame: nuclei of one frame share imaging conditions and the
#     annotators' selection, so they are not independent (28 frames -> enough levels).
both <- subset(w, matched_pre & matched_ft)
logit <- function(p) { p <- pmin(pmax(p, 1e-4), 1 - 1e-4); log(p / (1 - p)) }
both$diff_logit <- logit(both$iou_ft) - logit(both$iou_pre)
both$diff_raw   <- both$iou_ft - both$iou_pre

cat("\n=== Outline quality (nuclei detected by both models) ===\n")
cat(sprintf("n = %d nuclei in %d frames | mean raw IoU change %+.3f | improved in %.0f %% of nuclei\n",
            nrow(both), length(unique(both$frame)), mean(both$diff_raw), 100 * mean(both$diff_raw > 0)))

m1 <- lmer(diff_logit ~ 1 + (1 | frame), data = both, REML = TRUE)
print(summary(m1), correlation = FALSE)

ci_of <- function(model, parm) {
  # profile CI; falls back to Wald if profiling fails (e.g. singular fit)
  tryCatch(confint(model, parm = parm, method = "profile", quiet = TRUE),
           error = function(e) confint(model, parm = parm, method = "Wald"))
}
b1  <- fixef(m1)[["(Intercept)"]]
ci1 <- ci_of(m1, "(Intercept)")
vc  <- as.data.frame(VarCorr(m1))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\nMean change in logit(IoU): %+.3f (95%% CI %+.3f to %+.3f)\n", b1, ci1[1], ci1[2]))
cat(sprintf("Share of variance between frames (ICC): %.2f\n", icc))
cat(sprintf("Example: a nucleus with IoU 0.85 before -> %.3f after (CI %.3f-%.3f)\n",
            plogis(qlogis(0.85) + b1), plogis(qlogis(0.85) + ci1[1]), plogis(qlogis(0.85) + ci1[2])))

# Sanity check: ignoring frames (one-sample t-test). Similar result if frames add little.
tt <- t.test(both$diff_logit)
cat(sprintf("Check without frame effect (t-test): %+.3f (95%% CI %+.3f to %+.3f)\n",
            tt$estimate, tt$conf.int[1], tt$conf.int[2]))

# --- 4. Does the change differ near divisions? --------------------------------
b2 <- NA; ci2 <- c(NA, NA)
if (length(unique(both$near_division)) == 2) {
  m2 <- lmer(diff_logit ~ near_division + (1 | frame), data = both, REML = TRUE)
  b2  <- fixef(m2)[["near_divisionTRUE"]]
  ci2 <- ci_of(m2, "near_divisionTRUE")
  cat(sprintf("\nNear division (n = %d) vs. not (n = %d): extra change in logit(IoU) %+.3f (95%% CI %+.3f to %+.3f)\n",
              sum(both$near_division), sum(!both$near_division), b2, ci2[1], ci2[2]))
} else {
  cat("\nNo comparison near vs. not near division: only one group among nuclei detected by both models.\n")
}

# --- 5. Diagnostics -------------------------------------------------------------
png("results/figures/r_diagnostics.png", width = 1600, height = 1200, res = 150)
par(mfrow = c(2, 2))
plot(fitted(m1), resid(m1), xlab = "fitted", ylab = "residual",
     main = "Residuals vs fitted (m1)"); abline(h = 0, lty = 2)
qqnorm(resid(m1), main = "Normal Q-Q: residuals"); qqline(resid(m1))
re <- ranef(m1)$frame[, 1]
qqnorm(re, main = "Normal Q-Q: frame effects"); qqline(re)
hist(both$diff_logit, breaks = 30, main = "Paired differences", xlab = "logit(IoU) fine-tuned - pretrained")
invisible(dev.off())

png("results/figures/r_iou_change_by_division.png", width = 1000, height = 800, res = 150)
boxplot(diff_raw ~ near_division, data = both, names = c("not near", "near division (<= 1 h)"),
        ylab = "IoU change (fine-tuned - pretrained)", main = "Outline quality change per nucleus")
abline(h = 0, lty = 2)
invisible(dev.off())

# --- 6. Summary table -------------------------------------------------------------
summ <- data.frame(
  quantity = c("nuclei", "lost", "gained", "mcnemar_p",
               "logit_change", "logit_change_lo", "logit_change_hi", "icc_frame",
               "near_div_extra", "near_div_extra_lo", "near_div_extra_hi"),
  value = c(nrow(w), lost, gained, if (lost + gained > 0) mc$p.value else NA,
            b1, ci1[1], ci1[2], icc, b2, ci2[1], ci2[2])
)
write.csv(summ, "results/r_paired_comparison_summary.csv", row.names = FALSE)
cat("\nSaved: results/r_paired_comparison_summary.csv, results/figures/r_*.png\n")
