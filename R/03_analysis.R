# =============================================================================
# 03_analysis.R - all statistical analyses reported in the article
# Uses A, cls, part, flagged (02_prepare.R). Creates: res (per side), age, asym, me, sens_out
# =============================================================================

## ---------------------------------------------------------------- participants (Table 3)
desc <- function(a) c(n = length(a), mean = mean(a), sd = sd(a), median = median(a),
                      q1 = unname(quantile(a, .25)), q3 = unname(quantile(a, .75)), min = min(a), max = max(a))
part$Group <- factor(paste(part$Ethnicity, ifelse(part$Sex == "Female", "females", "males")),
                     levels = c("Fulani females", "Hausa females", "Fulani males", "Hausa males"))
age <- list()
age$tab <- rbind(t(sapply(split(part$Age, part$Group), desc)),
                 "All females" = desc(part$Age[part$Sex == "Female"]), "All males" = desc(part$Age[part$Sex == "Male"]),
                 Total = desc(part$Age))
mw <- function(f, sub = rep(TRUE, nrow(part))) wilcox.test(f, data = part[sub, ], exact = FALSE)$p.value
age$p_shapiro <- shapiro.test(residuals(aov(Age ~ Sex * Ethnicity, data = part)))$p.value
age$p_sex <- mw(Age ~ Sex); age$p_eth <- mw(Age ~ Ethnicity)
age$p_eth_f <- mw(Age ~ Ethnicity, part$Sex == "Female"); age$p_eth_m <- mw(Age ~ Ethnicity, part$Sex == "Male")
age$med <- c(sapply(split(part$Age, part$Sex), median), sapply(split(part$Age, part$Ethnicity), median))

## ---------------------------------------------------------------- analyses of each side
res <- list()
for (SIDE in c("Right", "Left")) {
  k <- cls$Side == SIDE; cs <- droplevels(cls[k, ])
  shape <- gpagen(A[, , k], print.progress = FALSE)$coords; M <- mshape(shape)
  ps <- pca_scores(shape); ev <- ps$d; pct <- 100 * ev / sum(ev); K <- npc_for(ev); X <- ps$x[, 1:K]
  bstick <- rev(cumsum(1 / rev(seq_along(ev)))) / length(ev) * 100
  tc <- tangent_check(A[, , k])
  pc_shape <- function(j, sdv) M + matrix(ps$pca$rotation[, j] * sdv * sqrt(ev[j]), ncol = 2, byrow = TRUE)
  pc_disp <- sapply(1:2, function(j) 100 * sqrt(rowSums((pc_shape(j, 2) - pc_shape(j, -2))^2)))

  ## Procrustes ANOVA (type II SS, residual randomization) with and without age; Pillai's trace on the PCs
  gdf <- geomorph.data.frame(shape = shape, Age = cs$Age, Sex = cs$Sex, Ethnicity = cs$Ethnicity)
  at <- anova_tab(procD.lm(shape ~ Age + Sex * Ethnicity, data = gdf, iter = ITER, SS.type = "II", print.progress = FALSE, seed = SEED))
  at0 <- anova_tab(procD.lm(shape ~ Sex * Ethnicity, data = gdf, iter = ITER, SS.type = "II", print.progress = FALSE, seed = SEED))
  man <- summary(manova(X ~ cs$Age + cs$Sex * cs$Ethnicity), test = "Pillai")$stats
  pill <- man[c("cs$Age", "cs$Sex", "cs$Ethnicity", "cs$Sex:cs$Ethnicity"), c("Pillai", "Pr(>F)")]
  ## simple effects (age as covariate)
  se <- do.call(rbind, lapply(list(c("Sex", "Ethnicity", "Fulani"), c("Sex", "Ethnicity", "Hausa"), c("Ethnicity", "Sex", "Female"),
                                   c("Ethnicity", "Sex", "Male")), function(z) {
    s <- cs[[z[2]]] == z[3]; g2 <- geomorph.data.frame(shape = shape[, , s], f = droplevels(cs[[z[1]]][s]), Age = cs$Age[s])
    t <- anova(procD.lm(shape ~ Age + f, data = g2, iter = ITER, SS.type = "II", print.progress = FALSE, seed = SEED))$table
    data.frame(Effect = z[1], Within = z[3], Rsq = 100 * t$Rsq[2], F = t$F[2], Z = t$Z[2], P = t$`Pr(>F)`[2]) }))

  ## CVA of the four groups with Morpho (pooled within-group covariance of all groups, 'ITER' permutations):
  ## Mahalanobis distances and classification on the PCs, Procrustes distances on all coordinates
  g4 <- cs$Group
  set.seed(SEED); mcv <- Morpho::CVA(X, g4, rounds = ITER, cv = TRUE, prior = rep(.25, 4), plot = FALSE)
  set.seed(SEED); mcvf <- Morpho::CVA(two.d.array(shape), g4, rounds = ITER, cv = FALSE, plot = FALSE)
  pr <- combn(levels(g4), 2); getd <- function(m, j) as.matrix(m)[pr[2, j], pr[1, j]]
  pw <- data.frame(Pair = paste(GLAB[pr[1, ]], "vs", GLAB[pr[2, ]]),
                   d = sapply(1:6, function(j) getd(mcvf$Dist$GroupdistEuclid, j)), P_d = sapply(1:6, function(j) getd(mcvf$Dist$probsEuclid, j)),
                   D = sapply(1:6, function(j) getd(mcv$Dist$GroupdistMaha, j)), P_D = sapply(1:6, function(j) getd(mcv$Dist$probsMaha, j)))
  t4 <- table(g4, lda(X, g4, prior = rep(.25, 4), CV = TRUE)$class)
  stopifnot("Morpho and MASS allocations differ" = identical(as.character(mcv$class), as.character(lda(X, g4, prior = rep(.25, 4), CV = TRUE)$class)))

  ## two-group DFA for sex and ethnicity; classification within subgroups; sensitivity to the number of PCs
  dsex <- two_group(X, cs$Sex); deth <- two_group(X, cs$Ethnicity)
  sub <- sapply(list(eth_f = c("Ethnicity", "Sex", "Female"), eth_m = c("Ethnicity", "Sex", "Male"),
                     sex_fu = c("Sex", "Ethnicity", "Fulani"), sex_ha = c("Sex", "Ethnicity", "Hausa")), function(z) {
    s <- cs[[z[2]]] == z[3]; pk <- gm.prcomp(shape[, , s]); Xs <- pk$x[, 1:npc_for(pk$d)]
    acc(table(cs[[z[1]]][s], lda(Xs, droplevels(cs[[z[1]]][s]), prior = c(.5, .5), CV = TRUE)$class)) })
  sens <- do.call(rbind, lapply(c(.90, .95, .99), function(pv) { kk <- npc_for(ev, pv); Xk <- ps$x[, 1:kk]
    f <- function(v) acc(table(v, lda(Xk, v, prior = rep(1 / nlevels(v), nlevels(v)), CV = TRUE)$class))
    data.frame(Variance = pv, PCs = kk, Sex = f(cs$Sex), Ethnicity = f(cs$Ethnicity), FourGroups = f(cs$Group)) }))

  ## differences between group means at each point (% of centroid size)
  mean_of <- function(f, l) mshape(shape[, , cs[[f]] == l])
  msex <- list(mean_of("Sex", "Female"), mean_of("Sex", "Male")); meth <- list(mean_of("Ethnicity", "Fulani"), mean_of("Ethnicity", "Hausa"))
  dd <- data.frame(Sex = 100 * sqrt(rowSums((msex[[2]] - msex[[1]])^2)), Ethnicity = 100 * sqrt(rowSums((meth[[2]] - meth[[1]])^2)))

  res[[SIDE]] <- list(shape = shape, cs = cs, M = M, ps = ps, pct = pct, bstick = bstick, K = K, X = X, tc = tc,
                      pc_shape = pc_shape, pc_disp = pc_disp, at = at, at0 = at0, pill = pill, se = se, pw = pw, t4 = t4,
                      dsex = dsex, deth = deth, sub = sub, sens = sens, dd = dd, msex = msex, meth = meth)
}

## ---------------------------------------------------------------- bilateral asymmetry (all participants)
asym <- list()
bs <- bilat.symmetry(A = A, ind = factor(cls$Individual), side = cls$Side, object.sym = FALSE, iter = ITER, seed = SEED, print.progress = FALSE)
st <- bs$shape.anova
asym$tab <- data.frame(Effect = rownames(st), df = st$Df, SS = st$SS, MS = st$MS, Rsq = 100 * st$Rsq, F = st$F, Z = st$Z, P = st$`Pr(>F)`)
Aall <- gpagen(A, print.progress = FALSE)$coords
idx <- split(seq_len(dim(A)[3]), cls$Individual)
asy <- simplify2array(lapply(idx, function(i) Aall[, , i[cls$Side[i] == "Right"]] - Aall[, , i[cls$Side[i] == "Left"]]))
ci <- cls[sapply(idx, function(i) i[cls$Side[i] == "Right"]), ]
asym$group <- anova_tab(procD.lm(asy ~ Age + Sex * Ethnicity, iter = ITER, SS.type = "II", print.progress = FALSE, seed = SEED,
                                 data = geomorph.data.frame(asy = asy, Age = ci$Age, Sex = ci$Sex, Ethnicity = ci$Ethnicity)))
asym$mR <- mshape(Aall[, , cls$Side == "Right"]); asym$mL <- mshape(Aall[, , cls$Side == "Left"]); asym$M <- mshape(Aall)
asym$pts <- 100 * sqrt(rowSums((asym$mR - asym$mL)^2)); asym$d <- sqrt(sum((asym$mR - asym$mL)^2))

## ---------------------------------------------------------------- measurement error (20 participants, both ears twice)
me <- list()
rep48 <- read48(RAW$replicates, renames); rid <- dimnames(rep48)[[3]]
stopifnot("re-digitized ear not in the main sample" = all(rid %in% dimnames(A)[[3]]))
isL <- substr(rid, 3, 3) == "L"; rep48[, 1, isL] <- -rep48[, 1, isL]
AM <- abind(A[, , rid], rep48, along = 3)
info <- cbind(classify_ids(c(rid, rid)), Replicate = factor(rep(1:2, each = length(rid))))
info$Individual <- factor(info$Individual)
SH <- gpagen(AM, curves = sliders, ProcD = FALSE, print.progress = FALSE)$coords
me$n_ind <- nlevels(info$Individual); me$per_group <- unique(table(substr(levels(info$Individual), 1, 2)))
for (sd in c("Right", "Left")) {
  k <- info$Side == sd; sh <- gpagen(SH[, , k], print.progress = FALSE)$coords; inf <- droplevels(info[k, ])
  gdf <- geomorph.data.frame(sh = sh, ind = inf$Individual, rp = inf$Replicate)
  t <- anova(procD.lm(sh ~ ind, data = gdf, iter = ITER, print.progress = FALSE, seed = SEED))$table
  ts <- anova(procD.lm(sh ~ ind + rp, data = gdf, iter = ITER, print.progress = FALSE, seed = SEED))$table
  MSi <- t$MS[1]; MSe <- t$MS[2]; s2i <- (MSi - MSe) / 2
  o1 <- which(inf$Replicate == 1); o2 <- which(inf$Replicate == 2)[match(inf$ID[o1], inf$ID[inf$Replicate == 2])]
  shift <- sapply(seq_along(o1), function(i) 100 * sqrt(rowSums((sh[, , o1[i]] - sh[, , o2[i]])^2)))
  me[[sd]] <- list(F = t$F[1], df1 = t$Df[1], df2 = t$Df[2], P = t$`Pr(>F)`[1], MSi = MSi, MSe = MSe,
                   err = 100 * t$SS[2] / sum(t$SS[1:2]), rep = s2i / (s2i + MSe),
                   sess = 100 * ts$Rsq[2], sess_P = ts$`Pr(>F)`[2], shift = apply(shift, 1, median))
}
grp <- c(paste0(1:11, " ", LM_NAMES), paste0("Curve C", rep(1:7, SEMI)))
me$pts <- aggregate(cbind(Right = me$Right$shift, Left = me$Left$shift), list(Element = factor(grp, levels = unique(grp))), mean)
bsm <- bilat.symmetry(A = SH, ind = info$Individual, side = info$Side, replicate = info$Replicate, object.sym = FALSE,
                      iter = ITER, seed = SEED, print.progress = FALSE)
sm <- bsm$shape.anova
me$asym <- data.frame(Effect = rownames(sm), df = sm$Df, SS = sm$SS, MS = sm$MS, Rsq = 100 * sm$Rsq, F = sm$F, Z = sm$Z, P = sm$`Pr(>F)`)

## ---------------------------------------------------------------- sensitivity: without ears flagged as shape outliers
sens_out <- lapply(c(Right = "Right", Left = "Left"), function(SIDE) {
  fl <- flagged[substr(flagged, 3, 3) == substr(SIDE, 1, 1)]
  if (!length(fl)) return(NULL)
  k <- cls$Side == SIDE & !cls$ID %in% fl
  sh <- gpagen(A[, , k], print.progress = FALSE)$coords; pk <- gm.prcomp(sh); Xs <- pk$x[, 1:npc_for(pk$d)]
  f <- function(v) acc(table(v, lda(Xs, v, prior = c(.5, .5), CV = TRUE)$class))
  t2 <- anova(procD.lm(sh ~ Age + Sex * Ethnicity, iter = 999, SS.type = "II", print.progress = FALSE, seed = SEED,
                       data = geomorph.data.frame(sh = sh, Age = cls$Age[k], Sex = cls$Sex[k], Ethnicity = cls$Ethnicity[k])))$table
  list(ids = fl, sex = f(droplevels(cls$Sex[k])), eth = f(droplevels(cls$Ethnicity[k])), r2_sex = 100 * t2$Rsq[2], r2_eth = 100 * t2$Rsq[3])
})
