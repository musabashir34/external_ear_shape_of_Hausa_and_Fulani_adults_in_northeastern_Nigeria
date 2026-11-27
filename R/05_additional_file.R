# =============================================================================
# 05_additional_file.R - Additional file 1 (Tables S1-S10): Word document, Excel workbook and CSV files
# =============================================================================
S <- list()
nm <- c(paste0(1:11, " ", LM_NAMES), paste0(12:48, " semilandmark (curve ", rep(1:7, SEMI), ")"))
rd <- function(x, d = 2) { x[] <- lapply(x, function(v) if (is.numeric(v)) round(v, d) else v); x }
S$S1 <- list("Median displacement between the two digitizations (% of centroid size), by landmark and curve (curve = mean of its semilandmarks)", rd(me$pts))
S$S2 <- list("Variance explained by the first ten principal components and broken-stick expectation (%)",
  rd(data.frame(PC = 1:10, Right_pct = res$Right$pct[1:10], Right_cum = cumsum(res$Right$pct)[1:10], Left_pct = res$Left$pct[1:10],
                Left_cum = cumsum(res$Left$pct)[1:10], BrokenStick = res$Right$bstick[1:10])))
S$S3 <- list("Displacement of each point between -2 SD and +2 SD along PC1 and PC2 (% of centroid size)",
  rd(data.frame(Point = nm, Right_PC1 = res$Right$pc_disp[, 1], Right_PC2 = res$Right$pc_disp[, 2], Left_PC1 = res$Left$pc_disp[, 1], Left_PC2 = res$Left$pc_disp[, 2]), 3))
S$S4 <- list("Correlation between Procrustes (geodesic) distances and Euclidean distances in tangent space for all pairs of ears",
  rd(rbind(cbind(Side = "Right", res$Right$tc), cbind(Side = "Left", res$Left$tc)), 5))
S$S5 <- list(sprintf("Procrustes ANOVA of ear shape without age (type II SS, %d permutations)", ITER),
  rd(rbind(cbind(Side = "Right", res$Right$at0), cbind(Side = "Left", res$Left$at0)), 4))
S$S6 <- list(sprintf("Simple effects of sex within ethnic groups and of ethnicity within sexes (age as covariate, %d permutations)", ITER),
  rd(rbind(cbind(Side = "Right", res$Right$se), cbind(Side = "Left", res$Left$se)), 4))
S$S7 <- list("Displacement of each point between group means and between sides (% of centroid size)",
  rd(data.frame(Point = nm, Right_sex = res$Right$dd$Sex, Right_ethnicity = res$Right$dd$Ethnicity, Left_sex = res$Left$dd$Sex,
                Left_ethnicity = res$Left$dd$Ethnicity, Directional_asymmetry = asym$pts), 3))
S$S8 <- list("Cross-validated classification of the four groups (rows: true group; columns: predicted group)",
  rbind(cbind(Side = "Right", Group = GLAB[rownames(res$Right$t4)], as.data.frame.matrix(res$Right$t4)),
        cbind(Side = "Left", Group = GLAB[rownames(res$Left$t4)], as.data.frame.matrix(res$Left$t4))))
S$S9 <- list("Cross-validated classification accuracy (%) by the proportion of shape variance retained",
  rd(rbind(cbind(Side = "Right", res$Right$sens), cbind(Side = "Left", res$Left$sens)), 1))
so <- do.call(rbind, lapply(names(Filter(Negate(is.null), sens_out)), function(s) with(sens_out[[s]],
  data.frame(Side = s, Excluded = paste(ids, collapse = ", "), Sex_CV_pct = sex, Ethnicity_CV_pct = eth, Sex_R2_pct = r2_sex, Ethnicity_R2_pct = r2_eth))))
S$S10 <- list("Results without the ears flagged as shape outliers", rd(so, 2))

dir.create(file.path(OUT, "additional_file_1"), showWarnings = FALSE, recursive = TRUE)
for (n in names(S)) write.csv(S[[n]][[2]], file.path(OUT, "additional_file_1", paste0(n, ".csv")), row.names = FALSE)
wb <- openxlsx::createWorkbook(); hs <- openxlsx::createStyle(textDecoration = "bold")
for (n in names(S)) {
  openxlsx::addWorksheet(wb, n)
  openxlsx::writeData(wb, n, paste0("Table ", n, ". ", S[[n]][[1]]), startRow = 1); openxlsx::addStyle(wb, n, hs, rows = 1, cols = 1)
  openxlsx::writeData(wb, n, S[[n]][[2]], startRow = 3, headerStyle = hs)
}
openxlsx::saveWorkbook(wb, file.path(OUT, "Additional_file_1.xlsx"), overwrite = TRUE)

## ---- Word version of Additional file 1 -------------------------------------
u8 <- function(x) { x <- as.character(x); Encoding(x) <- "UTF-8"; x }
LAB <- c(Element = "Landmark or curve", Right = "Right ear", Left = "Left ear", PC = "PC",
         Right_pct = "Right, %", Right_cum = "Right, cumulative %", Left_pct = "Left, %", Left_cum = "Left, cumulative %",
         BrokenStick = "Broken stick, %", Point = "Point", Right_PC1 = "Right PC1", Right_PC2 = "Right PC2",
         Left_PC1 = "Left PC1", Left_PC2 = "Left PC2", Side = "Side", Pairs = "Pairs of ears", r = "r", Slope = "Slope",
         Max_rho = "Maximum Procrustes distance", Effect = "Effect", df = "df", SS = "SS", MS = "MS", Rsq = "R² (%)",
         F = "F", Z = "Z", P = "P", Within = "Within", Right_sex = "Right, sex", Right_ethnicity = "Right, ethnicity",
         Left_sex = "Left, sex", Left_ethnicity = "Left, ethnicity", Directional_asymmetry = "Directional asymmetry",
         Group = "True group", FF = "Fulani female", FH = "Hausa female", MF = "Fulani male", MH = "Hausa male",
         Variance = "Shape variance retained (%)", PCs = "PCs", Sex = "Sex", Ethnicity = "Ethnicity", FourGroups = "Four groups",
         Excluded = "Ears excluded", Sex_CV_pct = "Sex, CV accuracy (%)", Ethnicity_CV_pct = "Ethnicity, CV accuracy (%)",
         Sex_R2_pct = "Sex, R² (%)", Ethnicity_R2_pct = "Ethnicity, R² (%)")
NOTE <- c(S1 = "Measured on the 20 participants digitized twice. Curve values are the mean over the semilandmarks of the curve.",
          S2 = "PC, principal component. Broken stick, variance expected for each component if variance were divided at random.",
          S3 = "Points 1–11 are fixed landmarks; points 12–48 are semilandmarks, numbered as in Table 2 of the article.",
          S4 = "r, Pearson correlation; slope of tangent-space on Procrustes distances; maximum Procrustes distance between any two ears (radians).",
          S5 = "Type II sums of squares. SS, sum of squares; MS, mean square; R², percentage of shape variation; Z, standardized effect size.",
          S6 = "Each row is a separate Procrustes ANOVA in the subgroup named under Within, with age as covariate. R², percentage of shape variation; Z, standardized effect size.",
          S7 = "Distance between the mean configurations of males and females, of Fulani and Hausa, and of right and reflected left ears.",
          S8 = "Discriminant functions on the principal components explaining 95% of shape variance; equal prior probabilities; leave-one-out cross-validation.",
          S9 = "Leave-one-out cross-validated accuracy of discriminant functions on the principal components (PCs) explaining the stated proportion of shape variance.",
          S10 = "Right ears whose Procrustes distance from the mean shape exceeded the upper quartile plus three times the interquartile range. CV, cross-validated.")
docx_tab <- function(n) {
  d <- S[[n]][[2]]
  if ("P" %in% names(d)) d$P <- fp(d$P)
  if ("Variance" %in% names(d)) d$Variance <- round(100 * d$Variance)
  d[] <- lapply(d, function(v) { v <- if (is.numeric(v)) format(v, scientific = FALSE, trim = TRUE, big.mark = "") else as.character(v)
                                  v[is.na(v) | v == "NA"] <- ""; u8(v) })
  if ("Effect" %in% names(d)) d$Effect <- sub(":", " \u00D7 ", d$Effect)
  ft <- flextable::flextable(d)
  hl <- vapply(names(d), function(k) if (k %in% names(LAB)) u8(LAB[[k]]) else k, "")
  ft <- flextable::set_header_labels(ft, values = as.list(hl))
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all"); ft <- flextable::fontsize(ft, size = 10, part = "all")
  ft <- flextable::bold(ft, part = "header"); ft <- flextable::theme_booktabs(ft); ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, j = seq_along(d)[-1], align = "center", part = "all")
  ft <- flextable::padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
  w1 <- if (names(d)[1] %in% c("Point", "Element")) 2.3 else if (names(d)[2] %in% c("Effect", "Group", "Excluded")) 0.6 else 0.8
  wr <- (6.5 - w1) / (ncol(d) - 1); if (names(d)[2] %in% c("Effect", "Group", "Excluded")) wr <- c(1.6, rep((6.5 - w1 - 1.6) / (ncol(d) - 2), ncol(d) - 2))
  ft <- flextable::width(ft, width = c(w1, rep_len(wr, ncol(d) - 1)))
  flextable::set_table_properties(ft, layout = "fixed")
}
tnr <- function(size = 12, bold = FALSE, italic = FALSE) officer::fp_text(font.family = "Times New Roman", font.size = size, bold = bold, italic = italic)
par1 <- function(..., align = "left") officer::fpar(..., fp_p = officer::fp_par(text.align = align, padding.bottom = 6))
TITLE_FALLBACK <- "Sexual dimorphism, ethnic variation and bilateral asymmetry in external ear shape of Hausa and Fulani adults in northeastern Nigeria"
doc <- officer::read_docx()
doc <- officer::body_add_fpar(doc, par1(officer::ftext("Supplementary material", tnr(14, TRUE))))
doc <- officer::body_add_fpar(doc, par1(officer::ftext(if (length(rmarkdown::metadata$title)) rmarkdown::metadata$title else TITLE_FALLBACK, tnr(12, italic = TRUE))))
doc <- officer::body_add_fpar(doc, par1(officer::ftext(paste0("Supplementary Tables S1–S", length(S),
         ". P values are from ", ITER, " permutations. Distances are given as a percentage of centroid size unless stated otherwise.") |> u8(), tnr(11))))
for (n in names(S)) {
  if (n != "S1") doc <- officer::body_add_break(doc)
  doc <- officer::body_add_fpar(doc, par1(officer::ftext(paste0("Table ", n, ". "), tnr(11, TRUE)), officer::ftext(u8(S[[n]][[1]]), tnr(11))))
  doc <- flextable::body_add_flextable(doc, docx_tab(n), align = "left")
  doc <- officer::body_add_fpar(doc, par1(officer::ftext(u8(NOTE[[n]]), tnr(9))))
}
print(doc, target = file.path(OUT, "Additional_file_1.docx"))
