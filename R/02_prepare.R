# =============================================================================
# 02_prepare.R - from the raw tpsDig2 files to the slid configurations
#   1. convert curves to 37 semilandmarks (48 points per ear), applying the image renamings
#   2. check the sample (400 right + 400 left ears, 100 per group, matching individuals)
#   3. reflect left ears and slide all 800 ears together (minimum bending energy)
#   4. screen for shape outliers (Procrustes distance from the side mean > Q3 + 3 IQR, after sliding)
# Creates: A (48 x 2 x 800 slid coordinates), cls (classifiers + age), flagged (outlier IDs)
# =============================================================================
renames <- if (file.exists(RAW$renames)) read.csv(RAW$renames, colClasses = "character") else NULL

if (file.exists(RAW$right) && file.exists(RAW$left)) {
  R48 <- read48(RAW$right, renames); L48 <- read48(RAW$left, renames)
  ids <- c(dimnames(R48)[[3]], dimnames(L48)[[3]])
  stopifnot("duplicated specimen IDs" = !anyDuplicated(ids), "right-ear file contains a left ear" = all(substr(dimnames(R48)[[3]], 3, 3) == "R"),
            "left-ear file contains a right ear" = all(substr(dimnames(L48)[[3]], 3, 3) == "L"))
  L48[, 1, ] <- -L48[, 1, ]                              # reflect left ears to right-ear orientation
  A0 <- abind(R48, L48, along = 3)
  A <- gpagen(A0, curves = sliders, ProcD = FALSE, print.progress = FALSE)$coords
  dimnames(A)[[3]] <- ids
  ## outlier screening: Procrustes distance from the mean shape of each side, after sliding
  screen <- function(X) { g <- gpagen(X, print.progress = FALSE)$coords; M <- mshape(g)
    d <- apply(g, 3, function(s) sqrt(sum((s - M)^2))); q <- quantile(d, c(.25, .75)); names(d)[d > q[2] + 3 * diff(q)] }
  side <- substr(ids, 3, 3)
  flagged <- c(screen(A[, , side == "R"]), screen(A[, , side == "L"]))
  dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
  write_tps(A, "data/derived/ears_slid.TPS"); writeLines(flagged, "data/derived/flagged_outliers.txt")
  DATA_SOURCE <- "raw tpsDig2 files"
} else {
  warning("Raw tpsDig2 files not found (", RAW$right, ", ", RAW$left, "); using the slid coordinates in data/derived/ears_slid.TPS")
  A <- readland.tps("data/derived/ears_slid.TPS", specID = "ID", warnmsg = FALSE)
  flagged <- readLines("data/derived/flagged_outliers.txt")
  DATA_SOURCE <- "data/derived/ears_slid.TPS"
}

cls <- classify_ids(dimnames(A)[[3]])
part <- read.csv(RAW$participants)
cls$Age <- part$Age[match(cls$Individual, part$Individual)]
stopifnot("age missing for some ears" = !anyNA(cls$Age),
          "each participant must have one right and one left ear" = all(table(cls$Individual, cls$Side) == 1),
          "100 participants per group expected" = all(table(cls$Group) == 200))
