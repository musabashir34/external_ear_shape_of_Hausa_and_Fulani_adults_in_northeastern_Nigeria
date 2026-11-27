# =============================================================================
# 04_figures.R - Figures 1-6, written to output/figures as
#   PDF (vector), TIFF (600 dpi, LZW; for submission) and PNG (300 dpi; for viewing)
# Widths follow Springer: 84 mm (single column) or 174 mm (double column).
# Ears are drawn upright in lateral view (aurale superior at the top, face to the right).
# Colours: Okabe-Ito palette (colour-blind safe); groups also differ by symbol and line type.
# Uses: res, asym (03_analysis.R), FIG (output folder), N_LM, chains, GLAB (01_functions.R)
# =============================================================================
FONT <- if (Sys.info()[["sysname"]] == "Linux") "Liberation Sans" else "Arial"   # Arial or its metric clone
GCOL <- c(FF = "#D55E00", FH = "#E69F00", MF = "#0072B2", MH = "#56B4E9")   # females warm, males cool
GPCH <- c(FF = 16, FH = 1, MF = 17, MH = 2)                                   # Fulani filled, Hausa open
GLTY <- c(FF = 1, FH = 2, MF = 1, MH = 2)                                     # Fulani solid, Hausa dashed
SEXCOL <- c(Female = "#D55E00", Male = "#0072B2"); ETHCOL <- c(Fulani = "#009E73", Hausa = "#CC79A7")
MAG <- 3                                                                      # magnification of mean differences
.alpha <- function(col, a) adjustcolor(col, alpha.f = a)
u8 <- function(x) { Encoding(x) <- "UTF-8"; x }                              # keep \u characters intact in any locale

## ---------------------------------------------------------------- output
save_fig <- function(name, w_mm, h_mm, draw) {
  dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
  w <- w_mm / 25.4; h <- h_mm / 25.4
  for (f in c("pdf", "tiff", "png")) {
    file <- file.path(FIG, paste0(name, ".", f))
    if (f == "pdf") grDevices::cairo_pdf(file, w, h, family = FONT)
    if (f == "tiff") tiff(file, w, h, units = "in", res = 600, compression = "lzw", type = "cairo", family = FONT)
    if (f == "png") png(file, w, h, units = "in", res = 300, type = "cairo", family = FONT)
    par(family = FONT, cex = 1, ps = 9, las = 1, tck = -0.015, mgp = c(1.9, 0.45, 0), lwd = 0.8)
    draw(); dev.off()
  }
}
panel <- function(lab, line = 0.4, cex = 1) mtext(lab, side = 3, adj = 0, line = line, font = 2, cex = cex)

## ---------------------------------------------------------------- orientation and drawing
## rotation (and a reflection for display only, if needed) that puts landmark 2 straight above landmark 4
## and the face (landmarks 1, 5, 10, 11) to the right; the same transformation is applied to every shape of a panel
ear_orientation <- function(M) {
  v <- M[2, ] - M[4, ]; th <- pi / 2 - atan2(v[2], v[1]); R <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  ctr <- colMeans(M); f <- function(X) sweep(X, 2, ctr) %*% t(R)
  Mr <- f(M); flip <- mean(Mr[c(1, 5, 10, 11), 1]) < Mr[3, 1]
  function(X) { Y <- f(X); if (flip) Y[, 1] <- -Y[, 1]; Y }
}
WIRE <- c(chains, list(c(9, 10, 11)))              # the seven curves plus the tragus line (9-10-11)
draw_ear <- function(X, col = "black", lwd = 1.3, lty = 1, pts = TRUE, fixed_cex = 0.7, semi_cex = 0.35, fill = "white") {
  for (ch in WIRE) lines(X[ch, ], col = col, lwd = lwd, lty = lty)
  if (pts) { points(X[-(1:N_LM), ], pch = 21, cex = semi_cex, col = col, bg = fill, lwd = 0.6)
             points(X[1:N_LM, ], pch = 21, cex = fixed_cex, col = col, bg = col, lwd = 0.6) }
}
label_landmarks <- function(X, col = "grey15", cex = 0.75, offset = 0.45) {
  pos <- c(4, 3, 2, 1, 4, 3, 2, 1, 4, 4, 3)       # 1 below, 2 left, 3 above, 4 right of each landmark
  text(X[1:N_LM, 1], X[1:N_LM, 2], labels = 1:N_LM, pos = pos, offset = offset, cex = cex, col = col, font = 2)
}
blank_ear <- function(XS, pad = 0.10) {
  r <- apply(do.call(rbind, XS), 2, range); e <- max(diff(r)) * pad
  plot(NA, xlim = r[, 1] + c(-1, 1) * e, ylim = r[, 2] + c(-1, 1) * e, asp = 1, axes = FALSE, xlab = "", ylab = "")
}
## thin-plate spline deformation grid from src to tgt
tps_fun <- function(src, tgt) {
  U <- function(r2) ifelse(r2 == 0, 0, r2 * log(r2))
  K <- U(as.matrix(dist(src))^2); P <- cbind(1, src)
  L <- rbind(cbind(K, P), cbind(t(P), matrix(0, 3, 3))); W <- solve(L, rbind(tgt, matrix(0, 3, 2)))
  function(Z) cbind(U(outer(rowSums(Z^2), rowSums(src^2), "+") - 2 * Z %*% t(src)), 1, Z) %*% W
}
draw_grid <- function(src, tgt, n = 18, col = "grey60", lwd = 0.45, pad = 0.10) {
  f <- tps_fun(src, tgt); r <- apply(src, 2, range); e <- diff(r) * pad
  xs <- seq(r[1, 1] - e[1], r[2, 1] + e[1], length = n); ys <- seq(r[1, 2] - e[2], r[2, 2] + e[2], by = diff(xs[1:2]))
  for (x in xs) lines(f(cbind(x, seq(min(ys), max(ys), length = 100))), col = col, lwd = lwd)
  for (y in ys) lines(f(cbind(seq(min(xs), max(xs), length = 100), y)), col = col, lwd = lwd)
}
ellipse95 <- function(x, y, level = 0.95, n = 200) {
  S <- cov(cbind(x, y)); m <- c(mean(x), mean(y)); r <- sqrt(qchisq(level, 2)); t <- seq(0, 2 * pi, length = n); e <- eigen(S)
  sweep(cbind(cos(t), sin(t)) %*% t(e$vectors %*% diag(sqrt(pmax(e$values, 0)))) * r, 2, m, "+")
}
roc_curve <- function(pos, score) {                # pos: TRUE for the second group; higher score = second group
  o <- order(score, decreasing = TRUE); y <- pos[o]; s <- score[o]
  tp <- cumsum(y) / sum(y); fp <- cumsum(!y) / sum(!y); keep <- c(diff(s) != 0, TRUE)
  list(fpr = c(0, fp[keep]), tpr = c(0, tp[keep]), auc = (sum(rank(score)[pos]) - sum(pos) * (sum(pos) + 1) / 2) / (sum(pos) * sum(!pos)))
}
SIDES <- c("Right", "Left")

## ---------------------------------------------------------------- Figure 1: configuration
M1 <- res$Right$M; o1 <- ear_orientation(M1); X1 <- o1(M1)
save_fig("Figure1_configuration", 84, 100, function() {
  par(mar = c(1.6, 0.3, 0.3, 0.3), xpd = NA)
  blank_ear(list(X1), 0.08)
  draw_ear(X1, col = "black", lwd = 1, fixed_cex = 1, semi_cex = 0.6)
  label_landmarks(X1, cex = 0.85, offset = 0.55)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  legend("bottom", horiz = TRUE, c("Fixed landmark", "Sliding semilandmark"), pch = 21, pt.bg = c("black", "white"),
         pt.cex = c(1, 0.7), bty = "n", cex = 0.9, text.width = c(0.36, 0.4))
})

## ---------------------------------------------------------------- Figure 2: PCA scatter (four groups)
save_fig("Figure2_PCA", 174, 92, function() {
  par(mfrow = c(1, 2), mar = c(3.2, 3.4, 1.6, 0.6), oma = c(0, 0, 3.4, 0))
  for (k in 1:2) { r <- res[[SIDES[k]]]; x <- r$ps$x[, 1]; y <- r$ps$x[, 2]; g <- as.character(r$cs$Group)
    els <- lapply(names(GCOL), function(s) ellipse95(x[g == s], y[g == s]))
    lim <- function(v, e) range(c(v, unlist(lapply(els, function(z) z[, e]))))
    plot(NA, xlim = lim(x, 1), ylim = lim(y, 2), asp = 1, axes = FALSE,
         xlab = sprintf("PC1 (%.1f%%)", r$pct[1]), ylab = sprintf("PC2 (%.1f%%)", r$pct[2]))
    abline(h = 0, v = 0, col = "grey85", lwd = 0.6)
    for (s in names(GCOL)) points(x[g == s], y[g == s], pch = GPCH[s], cex = 0.5, col = .alpha(GCOL[s], 0.5))
    for (j in seq_along(GCOL)) { s <- names(GCOL)[j]; lines(els[[j]], col = GCOL[s], lwd = 1.5, lty = GLTY[s])
      mx <- mean(x[g == s]); my <- mean(y[g == s])
      points(mx, my, pch = 3, cex = 1.3, lwd = 2.4, col = "white"); points(mx, my, pch = 3, cex = 1.2, lwd = 1.5, col = GCOL[s]) }
    axis(1, lwd = 0.8); axis(2, lwd = 0.8); box(lwd = 0.8); panel(sprintf("(%s) %s ear", letters[k], SIDES[k])) }
  n <- table(res$Right$cs$Group)[names(GCOL)]
  par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("top", ncol = 2, sprintf("%s (n = %d)", GLAB[names(GCOL)], n), pch = GPCH, lty = GLTY, col = GCOL, lwd = 1.5,
         pt.cex = 0.8, seg.len = 2.6, bty = "n", cex = 0.85, text.width = 0.25, x.intersp = 0.8, y.intersp = 1.1)
  legend("top", inset = c(0, 0.105), "+  group mean;  ellipses enclose 95% of each group", bty = "n", cex = 0.8, text.col = "grey30")
})

## ---------------------------------------------------------------- Figure 3: shape changes along PC1 and PC2
save_fig("Figure3_PC_shapes", 174, 128, function() {
  par(mfrow = c(2, 4), mar = c(0.3, 0.3, 1.8, 0.3), oma = c(0, 1.6, 0, 0)); i <- 0
  for (side in SIDES) { r <- res[[side]]; o <- ear_orientation(r$M); X0 <- o(r$M)
    ext <- lapply(1:2, function(j) lapply(c(-2, 2), function(s) o(r$pc_shape(j, s))))
    for (j in 1:2) for (k in 1:2) { i <- i + 1
      blank_ear(c(list(X0), unlist(ext, recursive = FALSE)), 0.05)
      draw_ear(X0, col = "grey65", lwd = 1, pts = FALSE)
      draw_ear(ext[[j]][[k]], col = "black", lwd = 1.1, fixed_cex = 0.55, semi_cex = 0.28)
      panel(sprintf("(%s) PC%d, %s2 SD", letters[i], j, u8(c("\u2212", "+"))[k]), cex = 0.95) }
    mtext(paste(side, "ear"), side = 2, outer = TRUE, at = c(0.75, 0.25)[match(side, SIDES)], line = 0.3, font = 2, las = 0) }
})

## ---------------------------------------------------------------- Figure 4: sex and ethnicity mean shapes with TPS grids
save_fig("Figure4_shape_differences", 174, 150, function() {
  layout(rbind(matrix(1:8, 2, byrow = TRUE), c(9, 9, 10, 10)), heights = c(1, 1, 0.13))
  par(mar = c(0.3, 0.3, 1.8, 0.3), oma = c(0, 1.6, 0, 0)); i <- 0
  for (side in SIDES) { r <- res[[side]]; o <- ear_orientation(r$M)
    for (fac in c("msex", "meth")) { cl <- if (fac == "msex") SEXCOL else ETHCOL
      A0 <- o(r[[fac]][[1]]); B0 <- o(r[[fac]][[2]]); Bm <- A0 + MAG * (B0 - A0)
      i <- i + 1; blank_ear(list(A0, Bm), 0.08)
      draw_ear(A0, col = cl[1], lwd = 1.2, fixed_cex = 0.6, semi_cex = 0.3); draw_ear(Bm, col = cl[2], lwd = 1.2, fixed_cex = 0.6, semi_cex = 0.3)
      if (i == 1) label_landmarks(A0, cex = 0.7, offset = 0.35)
      panel(sprintf("(%s) %s", letters[i], ifelse(fac == "msex", "Sex", "Ethnicity")), cex = 0.95)
      i <- i + 1; blank_ear(list(A0, Bm), 0.08); draw_grid(A0, Bm, n = 16)
      draw_ear(Bm, col = "black", lwd = 1.1, fixed_cex = 0.5, semi_cex = 0.25)
      panel(sprintf("(%s) %s to %s", letters[i], names(cl)[1], ifelse(fac == "msex", "male", "Hausa")), cex = 0.95) }
    mtext(paste(side, "ear"), side = 2, outer = TRUE, at = c(0.79, 0.33)[match(side, SIDES)], line = 0.3, font = 2, las = 0) }
  for (cl in list(SEXCOL, ETHCOL)) { par(mar = c(0, 0, 0, 0)); plot.new()
    legend("center", horiz = TRUE, c(paste(names(cl)[1], "mean"), u8(sprintf("%s mean (difference \u00D7%d)", names(cl)[2], MAG))),
           col = cl, lwd = 1.5, pch = 21, pt.bg = cl, pt.cex = 0.7, bty = "n", cex = 0.85) }
})

## ---------------------------------------------------------------- Figure 5: ROC curves of the cross-validated discriminant functions
ROC <- lapply(setNames(SIDES, SIDES), function(side) lapply(list(Sex = res[[side]]$dsex, Ethnicity = res[[side]]$deth),
                function(d) roc_curve(d$g == levels(d$g)[2], d$sc_cv)))
save_fig("Figure5_ROC", 174, 92, function() {
  par(mfrow = c(1, 2), mar = c(3.2, 3.4, 1.6, 0.8))
  cols <- c(Sex = "black", Ethnicity = "#009E73"); ltys <- c(Sex = 1, Ethnicity = 2)
  for (k in 1:2) { side <- SIDES[k]
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", asp = 1, axes = FALSE,
         xlab = u8("1 \u2212 specificity"), ylab = "Sensitivity")
    abline(0, 1, col = "grey75", lty = 2, lwd = 0.8)
    for (v in names(cols)) lines(ROC[[side]][[v]]$fpr, ROC[[side]][[v]]$tpr, col = cols[v], lty = ltys[v], lwd = 1.6, type = "s")
    axis(1, lwd = 0.8); axis(2, lwd = 0.8); box(lwd = 0.8); panel(sprintf("(%s) %s ear", letters[k], side))
    legend("bottomright", sprintf("%s (AUC = %.2f)", c("Sex", "Ethnicity"), sapply(ROC[[side]], `[[`, "auc")),
           col = cols, lty = ltys, lwd = 1.6, bty = "n", cex = 0.85, inset = 0.02) }
})

## ---------------------------------------------------------------- Figure 6: directional asymmetry
o6 <- ear_orientation(asym$M); L0 <- o6(asym$mL); R0 <- o6(asym$mR); Rm <- L0 + MAG * (R0 - L0)
ASYCOL <- c(Left = "#E69F00", Right = "black")
save_fig("Figure6_directional_asymmetry", 174, 92, function() {
  par(mfrow = c(1, 2), mar = c(2.2, 0.4, 1.8, 0.4))
  blank_ear(list(L0, Rm), 0.08)
  draw_ear(L0, col = ASYCOL["Left"], lwd = 1.3); draw_ear(Rm, col = ASYCOL["Right"], lwd = 1.3)
  label_landmarks(Rm, offset = 0.6); panel("(a)")
  blank_ear(list(L0, Rm), 0.08); draw_grid(L0, Rm); draw_ear(Rm, col = "black", lwd = 1.2, fixed_cex = 0.55, semi_cex = 0.3)
  panel("(b) Left (reflected) to right")
  par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  legend("bottom", horiz = TRUE, c("Left ear mean (reflected)", u8(sprintf("Right ear mean (difference \u00D7%d)", MAG))),
         col = ASYCOL, lwd = 1.5, pch = 21, pt.bg = ASYCOL, pt.cex = 0.7, bty = "n", cex = 0.85)
})
