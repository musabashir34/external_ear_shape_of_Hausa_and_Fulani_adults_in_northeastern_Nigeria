# =============================================================================
# 01_functions.R - landmark configuration, conversion of raw tpsDig2 files,
# statistical helpers and plotting helpers used by manuscript.Rmd
# =============================================================================

## -------------configuration---------------------
N_LM   <- 11
CURVES <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(6, 7), c(7, 8), c(8, 9))   # landmark pairs bounding each curve
SEMI   <- c(6, 8, 10, 3, 5, 4, 1)                                               # semilandmarks per curve (37)
LM_NAMES <- c("Otobasion superius", "Aurale superior", "Aurale posterior", "Subaurale", "Otobasion inferius",
              "End of inferior crus", "Antihelix concavity", "Antitragus", "Intertragic notch", "Tragus", "Incisura")
chains <- local({ idx <- N_LM; lapply(seq_along(CURVES), function(k) {
  s <- idx + seq_len(SEMI[k]); idx <<- idx + SEMI[k]; c(CURVES[[k]][1], s, CURVES[[k]][2]) }) })
## sliders for gpagen: each semilandmark slides between its two neighbours on the same curve
sliders <- do.call(rbind, lapply(chains, function(ch) t(sapply(2:(length(ch) - 1), function(j) ch[(j - 1):(j + 1)]))))
colnames(sliders) <- c("before", "slide", "after")
GLAB <- c(FF = "Fulani female", FH = "Hausa female", MF = "Fulani male", MH = "Hausa male")

## -------------raw tpsDig2 files -> 48-point configurations----------------
## raw image name -> specimen ID, after optional renaming (data/raw/image_renames.csv: raw_image, id)
image_id <- function(img, renames = NULL) {
  b <- basename(gsub("\\\\", "/", img))
  if (!is.null(renames) && b %in% renames$raw_image) return(renames$id[match(b, renames$raw_image)])
  gsub("[[:space:]]+", "_", trimws(sub("\\.[A-Za-z]+$", "", b)))
}
read_tps_curves <- function(infile, renames = NULL) {
  x <- readLines(infile, warn = FALSE, encoding = "latin1"); x <- x[trimws(x) != ""]
  starts <- grep("^LM=", trimws(x), ignore.case = TRUE); ends <- c(starts[-1] - 1, length(x))
  num <- function(s) as.numeric(strsplit(trimws(s), "[[:space:]]+")[[1]][1:2])
  lapply(seq_along(starts), function(r) {
    rec <- x[starts[r]:ends[r]]; trec <- trimws(rec); nlm <- as.integer(sub("LM=", "", trec[1], ignore.case = TRUE))
    lm <- matrix(t(sapply(rec[2:(1 + nlm)], num, USE.NAMES = FALSE)), ncol = 2)
    pts_i <- grep("^POINTS=", trec, ignore.case = TRUE)
    curves <- lapply(pts_i, function(i) { n <- as.integer(sub("POINTS=", "", trec[i], ignore.case = TRUE))
      matrix(t(sapply(rec[(i + 1):(i + n)], num, USE.NAMES = FALSE)), ncol = 2) })
    img <- sub("^IMAGE=", "", grep("^IMAGE=", rec, value = TRUE, ignore.case = TRUE)[1])   # untrimmed: keeps exact raw name
    list(lm = lm, curves = curves, id = image_id(img, renames), image = basename(gsub("\\\\", "/", img)))
  })
}
resample_inner <- function(P, k) {          # k equally spaced points along a polyline, excluding its ends
  seg <- sqrt(rowSums(diff(P)^2)); cum <- c(0, cumsum(seg)); at <- cum[length(cum)] * seq_len(k) / (k + 1)
  t(sapply(at, function(s) { j <- min(max(which(cum <= s)), nrow(P) - 1); f <- (s - cum[j]) / max(seg[j], 1e-9)
    P[j, ] + f * (P[j + 1, ] - P[j, ]) }))
}
## one ear: match each curve to its landmark pair, orient and anchor it, resample to the fixed number of semilandmarks
to48 <- function(r, max_gap = 0.08) {
  if (nrow(r$lm) != N_LM) stop(r$id, ": ", nrow(r$lm), " landmarks")
  size <- max(dist(r$lm)); cv <- r$curves[sapply(r$curves, nrow) >= 2]
  cost <- matrix(sapply(CURVES, function(p) sapply(cv, function(c) { A <- r$lm[p[1], ]; B <- r$lm[p[2], ]; s <- c[1, ]; e <- c[nrow(c), ]
    min(sqrt(sum((s - A)^2)) + sqrt(sum((e - B)^2)), sqrt(sum((s - B)^2)) + sqrt(sum((e - A)^2))) })), nrow = length(cv))
  out <- vector("list", length(CURVES)); used <- rep(FALSE, length(cv))
  for (step in seq_along(CURVES)) {
    m <- cost; m[used, ] <- Inf; m[, !sapply(out, is.null)] <- Inf; if (all(!is.finite(m))) break
    w <- which(m == min(m), arr.ind = TRUE)[1, ]; if (m[w[1], w[2]] / 2 > max_gap * size) break
    out[[w[2]]] <- cv[[w[1]]]; used[w[1]] <- TRUE
  }
  if (any(sapply(out, is.null))) stop(r$id, ": a curve could not be matched to its landmarks")
  semis <- lapply(seq_along(CURVES), function(k) { c <- out[[k]]; A <- r$lm[CURVES[[k]][1], ]; B <- r$lm[CURVES[[k]][2], ]
    if (sqrt(sum((c[1, ] - A)^2)) > sqrt(sum((c[1, ] - B)^2))) c <- c[nrow(c):1, , drop = FALSE]
    if (nrow(c) <= 2) stop(r$id, ": curve ", k, " has only its end points")
    resample_inner(rbind(A, c, B), SEMI[k]) })
  rbind(r$lm, do.call(rbind, semis))
}
read48 <- function(infile, renames = NULL) {  # raw tpsDig2 file -> 48 x 2 x n array, specimen IDs as names
  recs <- read_tps_curves(infile, renames); A <- simplify2array(lapply(recs, to48))
  dimnames(A)[[3]] <- sapply(recs, `[[`, "id"); A
}
write_tps <- function(A, file) {
  ids <- dimnames(A)[[3]]
  writeLines(unlist(lapply(seq_along(ids), function(i) c(paste0("LM=", dim(A)[1]), sprintf("%.8f %.8f", A[, 1, i], A[, 2, i]),
    paste0("IMAGE=", ids[i]), paste0("ID=", ids[i])))), file)
}
classify_ids <- function(ids) {             # "FFRE_004" -> sex, ethnicity, side, individual
  g <- substr(ids, 1, 2)
  data.frame(ID = ids, Side = factor(ifelse(substr(ids, 3, 3) == "R", "Right", "Left"), levels = c("Right", "Left")),
             Sex = factor(ifelse(substr(g, 1, 1) == "F", "Female", "Male")),
             Ethnicity = factor(ifelse(substr(g, 2, 2) == "F", "Fulani", "Hausa")),
             Group = factor(g, levels = c("FF", "FH", "MF", "MH")), Individual = paste0(g, "_", sub(".*_", "", ids)))
}

## ------------------ statistics--------------------------------
pca_scores <- function(shape) { p <- gm.prcomp(shape); keep <- p$d > max(p$d) * 1e-8; list(x = p$x[, keep], d = p$d[keep], pca = p) }
npc_for <- function(d, pv = PVAR) which(cumsum(d) / sum(d) >= pv)[1]
pooled_cov <- function(X, g) { g <- droplevels(g)
  r <- do.call(rbind, lapply(split(as.data.frame(X), g), function(z) scale(as.matrix(z), scale = FALSE)))
  crossprod(r) / (nrow(X) - nlevels(g)) }
mahal <- function(X, g, a, b, S) { d <- colMeans(X[g == a, , drop = FALSE]) - colMeans(X[g == b, , drop = FALSE]); sqrt(drop(t(d) %*% solve(S, d))) }
procd <- function(X, g, a, b) sqrt(sum((colMeans(X[g == a, , drop = FALSE]) - colMeans(X[g == b, , drop = FALSE]))^2))
perm_p <- function(obs, perm) (sum(perm >= obs) + 1) / (length(perm) + 1)
acc <- function(t) 100 * sum(diag(t)) / sum(t)
## tangent-space approximation: Procrustes (geodesic) distances vs Euclidean distances in tangent space
tangent_check <- function(Aside) {
  g1 <- gpagen(Aside, Proj = TRUE, print.progress = FALSE)$coords
  g0 <- gpagen(Aside, Proj = FALSE, print.progress = FALSE)$coords
  n <- dim(g0)[3]; dt <- as.vector(dist(t(apply(g1, 3, as.vector))))
  Z <- lapply(seq_len(n), function(i) { z <- scale(g0[, , i], scale = FALSE); z / sqrt(sum(z^2)) })
  rho <- unlist(lapply(1:(n - 1), function(i) vapply((i + 1):n, function(j)
    acos(min(1, sum(svd(crossprod(Z[[i]], Z[[j]]))$d))), numeric(1))))
  data.frame(Pairs = length(rho), r = cor(dt, rho), Slope = unname(coef(lm(dt ~ 0 + rho))), Max_rho = max(rho))
}
## two-group DFA: Procrustes distance and Mahalanobis D (permutation tests), Hotelling's T2, classification
two_group <- function(X, g, iter = ITER) {
  g <- droplevels(g); lv <- levels(g); n1 <- sum(g == lv[1]); n2 <- sum(g == lv[2]); p <- ncol(X)
  D <- mahal(X, g, lv[1], lv[2], pooled_cov(X, g)); dP <- procd(X, g, lv[1], lv[2])
  T2 <- (n1 * n2 / (n1 + n2)) * D^2; Fv <- T2 * (n1 + n2 - p - 1) / ((n1 + n2 - 2) * p)
  set.seed(SEED)
  perm <- replicate(iter, { gp <- sample(g); c(procd(X, gp, lv[1], lv[2]), mahal(X, gp, lv[1], lv[2], pooled_cov(X, gp))) })
  t_o <- table(g, predict(lda(X, g, prior = c(.5, .5)))$class); t_cv <- table(g, lda(X, g, prior = c(.5, .5), CV = TRUE)$class)
  sc_cv <- sapply(seq_len(nrow(X)), function(i) predict(lda(X[-i, ], g[-i], prior = c(.5, .5)), X[i, , drop = FALSE])$x[1])
  if (mean(sc_cv[g == lv[2]]) < 0) sc_cv <- -sc_cv
  list(d = dP, p_d = perm_p(dP, perm[1, ]), D = D, p_D = perm_p(D, perm[2, ]), T2 = T2, F = Fv, df1 = p, df2 = n1 + n2 - p - 1,
       p_T2 = pf(Fv, p, n1 + n2 - p - 1, lower.tail = FALSE), t_o = t_o, t_cv = t_cv, sc_cv = sc_cv, g = g)
}
anova_tab <- function(fit) { t <- anova(fit)$table
  data.frame(Effect = rownames(t), df = t$Df, SS = t$SS, MS = t$MS, Rsq = 100 * t$Rsq, F = t$F, Z = t$Z, P = t$`Pr(>F)`) }

## --------- number formatting for the text------------------
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
f4 <- function(x) formatC(x, format = "f", digits = 4)
## permutation P values: 4 decimals below 0.01, otherwise 3
fp <- function(p) ifelse(p < 0.01, formatC(p, format = "f", digits = 4), formatC(p, format = "f", digits = 3))
## parametric P values
fpp <- function(p) ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, format = "f", digits = 3)))
rng <- function(x, f = f1) paste0(f(min(x)), "\u2013", f(max(x)))

