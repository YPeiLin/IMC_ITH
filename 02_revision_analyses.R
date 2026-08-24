###############################################################################
##  Revision analysis
##
##  Ref1.1  robustness / independent validation
##            bootstrap + subsample stability of the Non-Cl Mo-ITH survival
##                association in the discovery cohort

##  Ref1.2  multivariable independence
##            Non-Cl Mo-ITH adjusted for clinicopathological variables, and
##            additionally for each established spatial metric

##  Ref1.4  redundancy of the prognostic metric vs established spatial metrics
##

###############################################################################

rm(list = ls())
suppressPackageStartupMessages(library(survival))
set.seed(20260821)

OUT   <- "sITH_comparison_out"
PRIM  <- "Non-Cl Mo_ITH"          # Select cell-type
AGG   <- "IMM2_ITH"               # the overall s-ITH, for contrast
NBOOT <- 2000
NSUB  <- 1000

log_lines <- character(0)
say <- function(...) { m <- paste0(...); cat(m, "\n"); log_lines <<- c(log_lines, m) }
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
g3  <- function(x) formatC(x, format = "g", digits = 3)

CMP <- c("Shannon_all", "Simpson_all", "Evenness_all", "Shannon_imm",
         "MixingScore", "MixingScore_norm", "TumWithImmNbr",
         "L_imm_r20", "L_imm_r50", "Lcross_immtum_r20", "Lcross_immtum_r50")
CLIN <- c("Sex_F", "Age_ge75", "BMI_ge30", "Smoker", "Stage_late")

## ---------------------------------------------------------------- load ----
cd <- readRDS(file.path(OUT, "cache_discovery.rds"))
cv <- readRDS(file.path(OUT, "cache_validation.rds"))
D  <- cbind(cd$ith, cd$base_m[rownames(cd$ith), , drop = FALSE])
V  <- cbind(cv$ith, cv$base_m[rownames(cv$ith), , drop = FALSE])
say("== Revision analyses ==")
say("discovery images: ", nrow(D), " | validation images: ", nrow(V))

clin <- read.csv("clin_discovery.csv", stringsAsFactors = FALSE)
clin$id <- sub("^LUAD_", "", clin$Key)
mm <- match(rownames(D), clin$id)
S <- data.frame(
  id    = rownames(D),
  time  = clin$Survival_time[mm],
  event = clin$Death[mm],
  prog  = clin$Progression[mm],
  clin[mm, CLIN],
  hist  = factor(clin$Histology[mm]),
  stringsAsFactors = FALSE)
for (m in c(PRIM, AGG, CMP)) S[[m]] <- D[[m]]

## z-score so hazard ratios are per SD and comparable across metrics
Z <- S
for (m in c(PRIM, AGG, CMP)) Z[[m]] <- as.numeric(scale(S[[m]]))
alias <- setNames(make.names(paste0("m_", c(PRIM, AGG, CMP))), c(PRIM, AGG, CMP))
for (m in names(alias)) Z[[alias[[m]]]] <- Z[[m]]
A <- function(x) unname(alias[x])

ok <- is.finite(Z$time) & Z$time > 0 & Z$event %in% c(0, 1) & is.finite(Z[[A(PRIM)]])
Z  <- Z[ok, ]
say("analysable for OS: n = ", nrow(Z), ", deaths = ", sum(Z$event))

cx <- function(rhs, dat = Z)
  tryCatch(coxph(as.formula(paste0("Surv(time, event) ~ ", rhs)), data = dat),
           error = function(e) NULL)
row1 <- function(fit, term) {
  s <- summary(fit)
  c(HR = s$conf.int[term, 1], lo = s$conf.int[term, 3],
    hi = s$conf.int[term, 4], p = s$coefficients[term, 5],
    C = unname(s$concordance[1]))
}

## ============ Ref1.2  multivariable independence ========================= ##
say("\n---- Ref1.2 : multivariable independence of ", PRIM, " ----")

models <- list(
  "unadjusted"                         = A(PRIM),
  "+ sex, age, BMI, smoking, stage"    = paste(c(A(PRIM), CLIN), collapse = " + "),
  "+ clinical + histology"             = paste(c(A(PRIM), CLIN, "hist"), collapse = " + "),
  "+ clinical + overall s-ITH"         = paste(c(A(PRIM), CLIN, A(AGG)), collapse = " + ")
)
mv <- do.call(rbind, lapply(names(models), function(nm) {
  f <- cx(models[[nm]]); if (is.null(f)) return(NULL)
  data.frame(model = nm, t(row1(f, A(PRIM))), n = f$n, events = f$nevent,
             stringsAsFactors = FALSE)
}))
for (i in seq_len(nrow(mv)))
  say(sprintf("   %-34s HR=%.2f [%.2f-%.2f] p=%-9.3g C=%.3f  (n=%d, d=%d)",
              mv$model[i], mv$HR[i], mv$lo[i], mv$hi[i], mv$p[i], mv$C[i],
              mv$n[i], mv$events[i]))

## additionally adjusted for each established spatial metric, one at a time
say("\n   adjusted for clinical variables AND each established spatial metric:")
mv2 <- do.call(rbind, lapply(CMP, function(m) {
  f <- cx(paste(c(A(PRIM), CLIN, A(m)), collapse = " + "))
  if (is.null(f)) return(NULL)
  r <- row1(f, A(PRIM)); rc <- row1(f, A(m))
  data.frame(comparator = m, prim_HR = r["HR"], prim_lo = r["lo"],
             prim_hi = r["hi"], prim_p = r["p"], cmp_HR = rc["HR"],
             cmp_p = rc["p"], C = r["C"], n = f$n, stringsAsFactors = FALSE)
}))
mv2$prim_q <- p.adjust(mv2$prim_p, "BH")
for (i in seq_len(nrow(mv2)))
  say(sprintf("   + %-20s %s HR=%.2f [%.2f-%.2f] p=%-9.3g (q=%.3g) | cmp p=%.3g",
              mv2$comparator[i], PRIM, mv2$prim_HR[i], mv2$prim_lo[i],
              mv2$prim_hi[i], mv2$prim_p[i], mv2$prim_q[i], mv2$cmp_p[i]))
say(sprintf("   summary: %s remains p<0.05 in %d/%d models; HR range %.2f-%.2f",
            PRIM, sum(mv2$prim_p < 0.05), nrow(mv2),
            min(mv2$prim_HR), max(mv2$prim_HR)))

## fully-loaded model: clinical + ALL established spatial metrics at once
f_all <- cx(paste(c(A(PRIM), CLIN, A(CMP)), collapse = " + "))
if (!is.null(f_all)) {
  r <- row1(f_all, A(PRIM))
  say(sprintf("\n   simultaneously adjusted for clinical + all %d spatial ",
              length(CMP)),
      sprintf("metrics: HR=%.2f [%.2f-%.2f] p=%.3g (n=%d)",
              r["HR"], r["lo"], r["hi"], r["p"], f_all$n))
}
write.csv(mv,  file.path(OUT, "R1_multivariable_models.csv"), row.names = FALSE)
write.csv(mv2, file.path(OUT, "R2_adjusted_for_each_spatial_metric.csv"),
          row.names = FALSE)

## ============ Ref1.4  redundancy of the prognostic metric ================ ##
say("\n---- Ref1.4 : redundancy of ", PRIM, " with established metrics ----")
red <- do.call(rbind, lapply(CMP, function(m) {
  cd_ <- suppressWarnings(cor.test(S[[PRIM]], S[[m]], method = "spearman"))
  vv  <- suppressWarnings(cor(V[[PRIM]], V[[m]], method = "spearman",
                             use = "pairwise.complete.obs"))
  data.frame(comparator = m, rho_disc = unname(cd_$estimate),
             p_disc = cd_$p.value, rho_val = vv, stringsAsFactors = FALSE)
}))
red$q_disc <- p.adjust(red$p_disc, "BH")
for (i in seq_len(nrow(red)))
  say(sprintf("   %-20s rho=%+.3f (q=%.3g) | validation rho=%+.3f",
              red$comparator[i], red$rho_disc[i], red$q_disc[i], red$rho_val[i]))
say(sprintf("   |rho| discovery: %.2f-%.2f (median %.2f); validation: %.2f-%.2f",
            min(abs(red$rho_disc)), max(abs(red$rho_disc)),
            median(abs(red$rho_disc)), min(abs(red$rho_val)),
            max(abs(red$rho_val))))
write.csv(red, file.path(OUT, "R3_redundancy_prognostic_metric.csv"),
          row.names = FALSE)

## univariable Cox for every comparator, for the side-by-side comparison
uni <- do.call(rbind, lapply(c(PRIM, AGG, CMP), function(m) {
  f <- cx(A(m)); if (is.null(f)) return(NULL)
  data.frame(metric = m, t(row1(f, A(m))), n = f$n, stringsAsFactors = FALSE)
}))
uni$q <- p.adjust(uni$p, "BH")
uni <- uni[order(uni$p), ]
write.csv(uni, file.path(OUT, "R4_univariable_all_metrics.csv"), row.names = FALSE)
say("\n   univariable Cox, sorted:")
for (i in seq_len(nrow(uni)))
  say(sprintf("   %-20s HR=%.2f [%.2f-%.2f] p=%-9.3g q=%-8.3g C=%.3f",
              uni$metric[i], uni$HR[i], uni$lo[i], uni$hi[i], uni$p[i],
              uni$q[i], uni$C[i]))

## nested LRT: does the prognostic metric add beyond each comparator, and vice versa
lrt <- do.call(rbind, lapply(CMP, function(m) {
  keep <- complete.cases(Z[, c(A(PRIM), A(m), CLIN)])
  Zc <- Z[keep, ]
  f_c <- cx(paste(c(A(m), CLIN), collapse = " + "), Zc)
  f_p <- cx(paste(c(A(PRIM), CLIN), collapse = " + "), Zc)
  f_b <- cx(paste(c(A(PRIM), A(m), CLIN), collapse = " + "), Zc)
  if (any(vapply(list(f_c, f_p, f_b), is.null, TRUE))) return(NULL)
  gp <- function(a) { x <- anova(a, f_b); x[[grep("^Pr", names(x))[1]]][2] }
  data.frame(comparator = m, LRT_p_add_prim = gp(f_c), LRT_p_add_cmp = gp(f_p),
             C_cmp = unname(summary(f_c)$concordance[1]),
             C_both = unname(summary(f_b)$concordance[1]),
             stringsAsFactors = FALSE)
}))
lrt$LRT_q_add_prim <- p.adjust(lrt$LRT_p_add_prim, "BH")
write.csv(lrt, file.path(OUT, "R5_nested_LRT.csv"), row.names = FALSE)
say("\n   nested LRT (all models include the clinical covariates):")
for (i in seq_len(nrow(lrt)))
  say(sprintf("   vs %-20s +%s p=%-9.3g (q=%.3g) | +cmp p=%-9.3g | C %.3f->%.3f",
              lrt$comparator[i], "prim", lrt$LRT_p_add_prim[i],
              lrt$LRT_q_add_prim[i], lrt$LRT_p_add_cmp[i], lrt$C_cmp[i],
              lrt$C_both[i]))
say(sprintf("   %s adds information beyond %d/%d comparators (q<0.05: %d)",
            PRIM, sum(lrt$LRT_p_add_prim < 0.05), nrow(lrt),
            sum(lrt$LRT_q_add_prim < 0.05)))

## ============ Ref1.1a  bootstrap and subsample stability ================= ##
say("\n---- Ref1.1 : robustness by resampling (", NBOOT, " bootstraps, ",
    NSUB, " 80% subsamples) ----")

boot_stat <- function(idx) {
  d <- Z[idx, ]
  f <- cx(A(PRIM), d)
  if (is.null(f)) return(c(NA, NA))
  s <- summary(f)
  c(s$conf.int[A(PRIM), 1], s$coefficients[A(PRIM), 5])
}
n <- nrow(Z)
bt <- t(vapply(seq_len(NBOOT), function(i) boot_stat(sample(n, n, replace = TRUE)),
               numeric(2)))
sb <- t(vapply(seq_len(NSUB), function(i)
  boot_stat(sample(n, floor(0.8 * n))), numeric(2)))
colnames(bt) <- colnames(sb) <- c("HR", "p")

rep_line <- function(M, label) {
  hr <- M[, "HR"]; pv <- M[, "p"]
  say(sprintf(paste0("   %-22s HR median %.2f, 95%%%% percentile interval ",
                     "%.2f-%.2f | HR>1 in %.1f%%%% | p<0.05 in %.1f%%%%"),
              label, median(hr, na.rm = TRUE),
              quantile(hr, 0.025, na.rm = TRUE), quantile(hr, 0.975, na.rm = TRUE),
              100 * mean(hr > 1, na.rm = TRUE), 100 * mean(pv < 0.05, na.rm = TRUE)))
  c(median = median(hr, na.rm = TRUE),
    lo = quantile(hr, 0.025, na.rm = TRUE)[[1]],
    hi = quantile(hr, 0.975, na.rm = TRUE)[[1]],
    pct_HR_gt1 = 100 * mean(hr > 1, na.rm = TRUE),
    pct_p_lt05 = 100 * mean(pv < 0.05, na.rm = TRUE))
}
rb <- rep_line(bt, "bootstrap")
rs <- rep_line(sb, paste0("80% subsample"))
write.csv(rbind(bootstrap = rb, subsample80 = rs),
          file.path(OUT, "R6_resampling_stability.csv"))

pdf(file.path(OUT, "R_F6_resampling_stability.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1))
for (z in list(list(bt, "Bootstrap (n resampled with replacement)"),
               list(sb, "80% subsample"))) {
  hr <- z[[1]][, "HR"]
  hist(hr, breaks = 40, col = "#B2182B55", border = "#B2182B",
       main = paste0(z[[2]], "\n", PRIM), xlab = "Hazard ratio per 1 SD",
       ylab = "resamples")
  abline(v = 1, col = "grey40", lty = 2, lwd = 2)
  abline(v = median(hr, na.rm = TRUE), col = "#2166AC", lwd = 2)
  legend("topright", c("HR = 1", "median"), col = c("grey40", "#2166AC"),
         lty = c(2, 1), lwd = 2, bty = "n", cex = 0.8)
}
invisible(dev.off())

## ============ Ref1.1b  independent validation cohort ===================== ##
say("\n---- Ref1.1 : independent multi-region validation cohort ----")
cvl <- read.csv("clin_validation.csv", stringsAsFactors = FALSE)
cvl$id <- sub("^LUAD_", "", cvl$Key)
Vd <- data.frame(id = rownames(V), V[, c(PRIM, AGG, CMP)], check.names = FALSE)
Vd$prog <- cvl$Progression[match(Vd$id, cvl$id)]
Vd$pid  <- sub("[AB]$", "", Vd$id)
say("   ", nrow(Vd), " images from ", length(unique(Vd$pid)), " patients; ",
    "no survival time is available for this cohort, so progression (a binary ",
    "outcome) is used.")

## patient-level metric = mean of the two regions (the manuscript's average s-ITH)
P <- aggregate(Vd[, c(PRIM, AGG)], by = list(pid = Vd$pid), FUN = mean,
               na.rm = TRUE)
P$prog <- Vd$prog[match(P$pid, Vd$pid)]
P <- P[is.finite(P[[PRIM]]) & P$prog %in% c(0, 1), ]
say("   patients analysable: ", nrow(P), " (", sum(P$prog), " progressed)")

lg <- glm(prog ~ scale(get(PRIM)), data = P, family = binomial)
sl <- summary(lg)$coefficients[2, ]
or <- exp(c(sl[1], sl[1] - 1.96 * sl[2], sl[1] + 1.96 * sl[2]))
say(sprintf("   %s vs progression: OR per 1 SD = %.2f [%.2f-%.2f], p = %.3g",
            PRIM, or[1], or[2], or[3], sl[4]))
wt <- wilcox.test(P[[PRIM]] ~ P$prog)
say(sprintf("   Wilcoxon progressed vs not: p = %.3g (medians %.3f vs %.3f)",
            wt$p.value, median(P[[PRIM]][P$prog == 1]),
            median(P[[PRIM]][P$prog == 0])))

## same test for progression in the discovery cohort, as a bridge
Zp <- Z[Z$prog %in% c(0, 1), ]
lgd <- glm(prog ~ get(A(PRIM)), data = Zp, family = binomial)
sd_ <- summary(lgd)$coefficients[2, ]
ord <- exp(c(sd_[1], sd_[1] - 1.96 * sd_[2], sd_[1] + 1.96 * sd_[2]))
say(sprintf("   discovery, same outcome (progression, n=%d, %d events): ",
            nrow(Zp), sum(Zp$prog)),
    sprintf("OR = %.2f [%.2f-%.2f], p = %.3g", ord[1], ord[2], ord[3], sd_[4]))

## ============ Ref1.1c  region A vs region B agreement ==================== ##
say("\n---- Ref1.1 : region A vs B agreement for the prognostic metric ----")
W <- reshape(Vd[, c("pid", "id", PRIM)], idvar = "pid",
             timevar = "id", direction = "wide")
ab <- do.call(rbind, lapply(split(Vd, Vd$pid), function(g) {
  if (nrow(g) != 2) return(NULL)
  g <- g[order(g$id), ]
  data.frame(pid = g$pid[1], A = g[[PRIM]][1], B = g[[PRIM]][2])
}))
ab <- ab[is.finite(ab$A) & is.finite(ab$B), ]
ct <- suppressWarnings(cor.test(ab$A, ab$B, method = "spearman"))
say(sprintf("   %d patients with both regions: Spearman rho = %.3f (p = %.3g)",
            nrow(ab), ct$estimate, ct$p.value))
say("   (low-to-moderate agreement is the manuscript's own point: a single ",
    "1 mm2 core does not fully represent the tumour, which is why s-ITH is ",
    "useful as a representativeness check.)")
write.csv(ab, file.path(OUT, "R7_regionAB_prognostic_metric.csv"),
          row.names = FALSE)

writeLines(log_lines, file.path(OUT, "R_revision_log.txt"))
say("\n== done ==")
