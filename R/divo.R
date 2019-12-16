# divo.R
# Version: 1.0.1
# Autor: Christoph Sadee, Maciej Pietrzak, Michal Seweryn, Cankun Wang, Grzegorz Rempala
# Maintainer: Maciej Pietrzak <pietrzak.20@osu.edu>
# License: GPL (>=3)



cvg <- function(x) {
  y <- x
  n <- sum(y)
  f1 <- sum(y == 1)
  if (f1 == n) {
    f1 <- n - 1
  }
  return(1 - f1 / n)
}

i.inp <- function(x, alpha = 1, CI = 0.95, resample = 100, graph = FALSE,
                  csv_output = FALSE, PlugIn = FALSE, size = 1, CVG = FALSE,
                  saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (alpha > 0) {
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    if (0 < CI & CI < 1) {
      dimX <- dim(xx)
      OutMean <- mat.or.vec(dimX[2], dimX[2])
      OutMin <- mat.or.vec(dimX[2], dimX[2])
      OutMax <- mat.or.vec(dimX[2], dimX[2])
      OutCvg <- mat.or.vec(dimX[2], 1)
      .Call(
        "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "INP", alpha,
        resample, CI, 0, PlugIn, size, 0, CVG, saveBootstrap, OutMean, OutMin,
        OutMax, OutCvg
      )
      return(.e.clust(xx, csv_output, graph,
        funct = "INP", CVG, PlugIn, OutMean,
        OutMin, OutMax, OutCvg, saveBootstrap
      ))
    }
    else {
      stop("Confidence interval must be between 0 and 1; default CI=0.95")
    }
  }
  else {
    stop("alpha must be greater than 0")
  }
}

li <- function(x, CI = 0.95, resample = 100, graph = FALSE, csv_output = FALSE,
               PlugIn = FALSE, size = 1, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (class(xx)[1] != "matrix") {
    xx <- as.matrix(xx)
  }
  if (0 < CI & CI < 1) {
    dimX <- dim(xx)
    OutMean <- mat.or.vec(dimX[2], dimX[2])
    OutMin <- mat.or.vec(dimX[2], dimX[2])
    OutMax <- mat.or.vec(dimX[2], dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "LI", 0, resample,
      CI, 0, PlugIn, size, 0, 0, saveBootstrap, OutMean, OutMin, OutMax, OutCvg
    )
    return(.e.clust(xx, csv_output, graph,
      funct = "LI", FALSE, PlugIn, OutMean,
      OutMin, OutMax, OutCvg, saveBootstrap
    ))
  }
  else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
}

ji <- function(x, CI = 0.95, resample = 100, graph = FALSE, csv_output = FALSE,
               PlugIn = FALSE, size = 1, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (class(xx)[1] != "matrix") {
    xx <- as.matrix(xx)
  }
  if (0 < CI & CI < 1) {
    dimX <- dim(xx)
    OutMean <- mat.or.vec(dimX[2], dimX[2])
    OutMin <- mat.or.vec(dimX[2], dimX[2])
    OutMax <- mat.or.vec(dimX[2], dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "JI", 0, resample,
      CI, 0, PlugIn, size, 0, 0, saveBootstrap, OutMean, OutMin, OutMax, OutCvg
    )
    return(.e.clust(xx, csv_output, graph,
      funct = "JI", FALSE, PlugIn, OutMean,
      OutMin, OutMax, OutCvg, saveBootstrap
    ))
  }
  else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
}

rd <- function(x, alpha = 0.5, CI = 0.95, resample = 100, graph = FALSE,
               csv_output = FALSE, PlugIn = FALSE, size = 1, CVG = FALSE,
               saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (class(xx)[1] != "matrix") {
    xx <- as.matrix(xx)
  }
  chck <- 0
  n <- ncol(xx)
  for (a in seq(1, n, 1)) {
    if (a < n) {
      for (b in seq(a + 1, n, 1)) {
        if (sum(x[, a] * x[, b]) == 0) {
          chck <- chck + 1
          warning(as.character(paste("Populations: ", as.character(colnames(x)[a]), ",
                               ", as.character(colnames(x)[b]),
            " are orthogonal.",
            sep = ""
          )))
        }
      }
    }
  }
  if (chck == 0) {
    if (alpha > 0) {
      if (0 < CI & CI < 1) {
        dimX <- dim(xx)
        OutMean <- mat.or.vec(dimX[2], dimX[2])
        OutMin <- mat.or.vec(dimX[2], dimX[2])
        OutMax <- mat.or.vec(dimX[2], dimX[2])
        OutCvg <- mat.or.vec(dimX[2], 1)
        .Call(
          "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "RD", alpha,
          resample, CI, 0, PlugIn, size, 0, CVG, saveBootstrap, OutMean, OutMin,
          OutMax, OutCvg
        )
        return(.e.clust(xx, csv_output, graph,
          funct = "RD", CVG, PlugIn, OutMean,
          OutMin, OutMax, OutCvg, saveBootstrap
        ))
      }
      else {
        stop("Confidence interval must be between 0 and 1; default CI=0.95")
      }
    }
    else {
      stop("alpha must be greater than 0")
    }
  }
  else {
    stop("alpha must be != 1")
  }
}

srd <- function(x, alpha = 0.5, CI = 0.95, resample = 100, graph = FALSE,
                csv_output = FALSE, PlugIn = FALSE, size = 1, CVG = FALSE,
                saveBootstrap = FALSE) {
  if (alpha < 1) {
    xx <- x
    resample <- round(resample, 0)
    xx[is.na(xx)] <- 0
    if (class(xx)[1] == "numeric") {
      stop("Number of columns must be greater than 1")
    }
    if (resample < 1) {
      stop("resample must be greater than 1")
    }
    if (size < 0.1) {
      stop("size must be greater than or equal to 0.1")
    }
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    chck <- 0
    n <- ncol(xx)
    for (a in seq(1, n, 1)) {
      if (a < n) {
        for (b in seq(a + 1, n, 1)) {
          if (sum(x[, a] * x[, b]) == 0) {
            chck <- chck + 1
            warning(as.character(paste("Populations: ", as.character(colnames(x)[a]),
              ", ", as.character(colnames(x)[b]),
              " are orthogonal.",
              sep = ""
            )))
          }
        }
      }
    }
    if (chck == 0) {
      if (alpha > 0) {
        if (0 < CI & CI < 1) {
          if (class(xx)[1] != "matrix") {
            xx <- as.matrix(xx)
          }
          dimX <- dim(xx)
          OutMean <- mat.or.vec(dimX[2], dimX[2])
          OutMin <- mat.or.vec(dimX[2], dimX[2])
          OutMax <- mat.or.vec(dimX[2], dimX[2])
          OutCvg <- mat.or.vec(dimX[2], 1)
          .Call(
            "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "RDS", alpha,
            resample, CI, 0, PlugIn, size, 0, CVG, saveBootstrap, OutMean, OutMin,
            OutMax, OutCvg
          )
          return(.e.clust(xx, csv_output, graph,
            funct = "RD", CVG, PlugIn, OutMean,
            OutMin, OutMax, OutCvg, saveBootstrap
          ))
        }
        else {
          stop("Confidence interval must be between 0 and 1; default CI=0.95")
        }
      }
      else {
        stop("alpha must be greater than 0")
      }
    }
    else {
      stop("alpha must be < 1")
    }
  }
}

mh <- function(x, CI = 0.95, resample = 100, graph = FALSE, csv_output = FALSE,
               PlugIn = FALSE, size = 1, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (0 < CI & CI < 1) {
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    dimX <- dim(xx)
    OutMean <- mat.or.vec(dimX[2], dimX[2])
    OutMin <- mat.or.vec(dimX[2], dimX[2])
    OutMax <- mat.or.vec(dimX[2], dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "MH", 0, resample,
      CI, 0, PlugIn, size, 0, 0, saveBootstrap, OutMean, OutMin, OutMax, OutCvg
    )
    return(.e.clust(xx, csv_output, graph,
      funct = "MH", FALSE, PlugIn,
      OutMean, OutMin, OutMax, OutCvg, saveBootstrap
    ))
  }
  else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
}

pg <- function(x, alpha = 1, beta = alpha, CI = 0.95, resample = 100, graph = FALSE,
               csv_output = FALSE, PlugIn = FALSE, size = 1, CVG = FALSE,
               saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (alpha > 0 & beta > 0) {
    if (0 < CI & CI < 1) {
      if (class(xx)[1] != "matrix") {
        xx <- as.matrix(xx)
      }
      dimX <- dim(xx)
      OutMean <- mat.or.vec(dimX[2], dimX[2])
      OutMin <- mat.or.vec(dimX[2], dimX[2])
      OutMax <- mat.or.vec(dimX[2], dimX[2])
      OutCvg <- mat.or.vec(dimX[2], 1)
      .Call(
        "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "PG", alpha,
        resample, CI, 0, PlugIn, size, beta, CVG, saveBootstrap, OutMean, OutMin,
        OutMax, OutCvg
      )
      return(.e.clust(xx, csv_output, graph,
        funct = "PG", CVG, PlugIn, OutMean,
        OutMin, OutMax, OutCvg, saveBootstrap
      ))
    }
    else {
      stop("Confidence interval must be between 0 and 1; default CI=0.95")
    }
  }
  else {
    stop("alpha and beta must be greater than 0")
  }
}

pg.ht <- function(x, alpha = 1, beta = alpha, CI = 0.95, resample = 100, graph = FALSE,
                  csv_output = FALSE, PlugIn = FALSE, size = 1, CVG = FALSE,
                  saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (alpha > 0 & beta > 0) {
    if (0 < CI & CI < 1) {
      if (class(xx)[1] != "matrix") {
        xx <- as.matrix(xx)
      }
      dimX <- dim(xx)
      OutMean <- mat.or.vec(dimX[2], dimX[2])
      OutMin <- mat.or.vec(dimX[2], dimX[2])
      OutMax <- mat.or.vec(dimX[2], dimX[2])
      OutCvg <- mat.or.vec(dimX[2], 1)
      .Call(
        "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "PG_HT", alpha,
        resample, CI, 0, PlugIn, size, 0, CVG, saveBootstrap, OutMean, OutMin,
        OutMax, OutCvg
      )
      # .Call("read",xx,"PG_HT",alpha,resample,CI,"fa",PlugIn,size,1,CVG,
      # saveBootstrap,OutMean,OutMin,OutMax, OutCvg)
      return(.e.clust(xx, csv_output, graph,
        funct = "PG.ht", CVG, PlugIn,
        OutMean, OutMin, OutMax, OutCvg, saveBootstrap
      ))
    }
    else {
      stop("Confidence interval must be between 0 and 1; default CI=0.95")
    }
  }
  else {
    stop("alpha and beta must be greater than 0")
  }
}

i.in <- function(x, alpha = 1, CI = 0.95, resample = 100, PlugIn = FALSE, size = 1,
                 CVG = FALSE, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (class(xx)[1] == "numeric") {
    stop("Number of columns must be greater than 1")
  }
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (0.1 <= alpha & alpha <= 1) {
    xx <- x
    if (0 < CI & CI < 1) {
      if (class(xx)[1] != "matrix") {
        xx <- as.matrix(xx)
      }
      dimX <- dim(xx)
      OutMean <- matrix(0)
      OutMin <- matrix(0)
      OutMax <- matrix(0)
      OutCvg <- mat.or.vec(dimX[2], 1)
      .Call(
        "read", "OL", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, xx, "IN", alpha,
        resample, CI, 0, PlugIn, size, 0, CVG, saveBootstrap, OutMean, OutMin,
        OutMax, OutCvg
      )
      get.data <- rbind(OutMean, OutMin, OutMax)
      colnames(get.data) <- c("I.index")
      rownames(get.data) <- c("Mean", "Lower.Quantile", "Upper.Quantile")
      return(get.data)
    }
    else {
      stop("Confidence interval must be between 0 and 1; default CI=0.95")
    }
  }
  else {
    stop("alpha must be between 0.1 and 1")
  }
}

.e.clust <- function(xy, csv_output, graph, funct, CVG, PlugIn, OutMean, OutMin,
                     OutMax, OutCvg, saveBootstrap) {
  pr_mean <- OutMean
  pr_min <- OutMin
  pr_max <- OutMax
  colnames(pr_mean) <- colnames(xy) -> rownames(pr_mean)
  colnames(pr_min) <- colnames(xy) -> rownames(pr_min)
  colnames(pr_max) <- colnames(xy) -> rownames(pr_max)
  try(if (funct != "RD") {
    diag(pr_mean) <- 1
    diag(pr_min) <- 1
    diag(pr_max) <- 1
  }, silent = TRUE)
  outlist <- list(pr_mean, pr_min, pr_max)
  names(outlist) <- c("Mean", "Lower.Quantile", "Upper.Quantile")
  if (graph != TRUE & graph != FALSE) {
    fname <- paste(graph, ".pdf", sep = "")
  }
  if (graph == TRUE) {
    fname <- "divo_Clust.pdf"
  }
  if (graph != FALSE) {
    .p.clust(xy, fname, funct,
      lab = xy[0, ], PlugIn, OutMean,
      OutMin, OutMax
    )
  }
  if (csv_output != FALSE & csv_output != TRUE) {
    fname.csv <- paste(csv_output, ".csv", sep = "")
  }
  if (csv_output == TRUE) {
    fname.csv <- "divo_Overlap_out.csv"
  }
  if (csv_output != FALSE) {
    f.list.to.df <- function(csv.x) {
      colnames(csv.x[[1]]) -> cnames
      colnames(csv.x[[1]]) <- NULL
      rownames(csv.x[[1]]) <- NULL
      temp <- rep(0, length(csv.x[[1]][1, ]))
      for (i in 1:length(csv.x)) {
        colnames(csv.x[[i]]) <- NULL
        rownames(csv.x[[i]]) <- NULL
        temp <- rbind(temp, rep("", length(csv.x[[1]][1, ])), csv.x[[i]])
      }
      temp <- temp[-1, ]
      if (PlugIn == TRUE) {
        temp <- cbind(c(
          "PlugIn", cnames, "Lower.Quantile", cnames,
          "Upper.Quantile", cnames
        ), temp)
      }
      else {
        temp <- cbind(c(
          "Mean", cnames, "Lower.Quantile", cnames, "Upper.Quantile",
          cnames
        ), temp)
      }
      colnames(temp) <- c("Type", cnames)
      temp[temp == 0] <- ""
      data.frame(temp) -> temp
      rownames(temp) <- NULL
      return(temp)
    }
    f.list.to.df(outlist) -> overlap.list.to.df
    if (PlugIn == TRUE) {
      overlap.list.to.df <- overlap.list.to.df[1:dim(overlap.list.to.df)[2], ]
    }
    write.csv(file = fname.csv, x = overlap.list.to.df, row.names = FALSE)
  }

  if (PlugIn == TRUE) {
    outlist <- outlist[1]
    names(outlist)[length(outlist)] <- "PlugIn"
  }
  if (CVG == TRUE) {
    a <- as.data.frame(OutCvg)
    colnames(a) <- c("CVG")
    rownames(a) <- colnames(xy)
    outlist[[length(outlist) + 1]] <- a
    names(outlist)[length(outlist)] <- "Coverage"
  }
  if (saveBootstrap != TRUE & saveBootstrap != FALSE) {
    file.rename("divoBootstrapResults.txt", paste(saveBootstrap, "txt",
      sep = "."
    ))
  }
  return(outlist)
}

.p.clust <- function(xx, fname, funct, lab, PlugIn, OutMean = OutMean,
                     OutMin = OutMin, OutMax = OutMax, h = FALSE, height = 7,
                     width = 0.6 * height) {
  if (funct != "RD") {
    mat.1 <- rbind(lab, 1 - as.matrix(round(OutMean, 6)))
    mat.2 <- rbind(lab, 1 - as.matrix(round(OutMin, 6)))
    mat.3 <- rbind(lab, 1 - as.matrix(round(OutMax, 6)))
  }
  else {
    mat.1 <- rbind(lab, as.matrix(round(OutMean, 6)))
    mat.2 <- rbind(lab, as.matrix(round(OutMin, 6)))
    mat.3 <- rbind(lab, as.matrix(round(OutMax, 6)))
  }
  as.dendrogram(agnes(mat.1, diss = TRUE, method = "ward")) -> dend.mat.1
  as.dendrogram(agnes(mat.2, diss = TRUE, method = "ward")) -> dend.mat.2
  as.dendrogram(agnes(mat.3, diss = TRUE, method = "ward")) -> dend.mat.3
  cex.val <- 0.4
  if (ncol(xx) < 5) {
    cex.val <- 0.5
  }
  if (ncol(xx) > 30) {
    cex.val <- 0.3
  }
  if (PlugIn == FALSE) {
    par(mfrow = c(3, 1), cex = cex.val, lwd = 0.5, mar = c(10, 4, 2, 0.5))
    plot(dend.mat.1, horiz = h, main = "Mean")
    plot(dend.mat.2, horiz = h, main = "Lower Quantile")
    plot(dend.mat.3, horiz = h, main = "Upper Quantile")
    dev.copy(pdf, fname, height = height, width = width)
    dev.off()
  }
  if (PlugIn == TRUE) {
    par(mfrow = c(1, 1), mar = c(10, 4, 2, 0.5))
    plot(dend.mat.1, horiz = h, main = "PlugIn")
    dev.copy(pdf, fname)
    dev.off()
  }
}

dp <- function(x, alpha = seq(0.1, 2, 0.1), CI = 0.95, resample = 100,
               single_graph = FALSE, pooled_graph = FALSE, csv_output = FALSE,
               PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (0 < CI & CI < 1) {
    f <- "RE"
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    Alpha.Profile <- alpha
    if (min(Alpha.Profile) <= 0) {
      stop("alpha must be larger than 0")
    }
    DP_nAP <- length(Alpha.Profile)
    dimX <- dim(xx)
    OutMean <- mat.or.vec(DP_nAP, dimX[2])
    OutMin <- mat.or.vec(DP_nAP, dimX[2])
    OutMax <- mat.or.vec(DP_nAP, dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "DP", xx, "DP", 1, resample, CI, "RE", PlugIn, size,
      Alpha.Profile, CVG, saveBootstrap, OutMean, OutMin, OutMax,
      OutCvg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    if (csv_output != FALSE) {
      .csv.save(Alpha.Profile, OutMean, OutMin, OutMax,
        PlugIn, csv_output, CVG,
        def.fname.csv = "DP_output.csv", xx
      )
    }
    if (PlugIn == FALSE) {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(
          alpha = Alpha.Profile,
          Mean = OutMean[, l],
          Lower.Quantile = OutMin[, l],
          Upper.Quantile = OutMax[, l]
        )
      }
      names(li) <- colnames(xx)
    }
    else {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(Alpha = Alpha.Profile, PlugIn = OutMean[, l])
      }
      names(li) <- colnames(xx)
    }
    if (length(Alpha.Profile) <= 1) {
      warning("For graphics alpha must be a vector of lenght greater than 1")
    }
    y.label <- "Diversity Index"
    title <- "Diversity Profile:"
    single.fname <- "DP_"
    max.y <- (max(OutMean) + .3)
    min.y <- min(OutMean) - .3
    if (pooled_graph != FALSE & length(alpha) > 1) {
      .pooled.graph(Alpha.Profile, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, pooled_graph, PlugIn,
        def.fname.pdf = "divo_DP_AllPopulations.pdf", xx
      )
    }
    if (single_graph != FALSE & length(alpha) > 1) {
      .single.graph(
        Alpha.Profile, xx, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, single_graph, PlugIn, single.fname
      )
    }
  } else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
  return(li)
}

ens <- function(x, alpha = seq(0.1, 2, 0.1), CI = 0.95, resample = 100,
                single_graph = FALSE, pooled_graph = FALSE, csv_output = FALSE,
                PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (0 < CI & CI < 1) {
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    Alpha.Profile <- alpha
    if (min(Alpha.Profile) <= 0) {
      stop("alpha must be larger than 0")
    }
    DP_nAP <- length(Alpha.Profile)
    dimX <- dim(xx)
    OutMean <- mat.or.vec(DP_nAP, dimX[2])
    OutMin <- mat.or.vec(DP_nAP, dimX[2])
    OutMax <- mat.or.vec(DP_nAP, dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "DP", xx, "ENS", 1, resample, CI, "RE", PlugIn, size,
      Alpha.Profile, CVG, saveBootstrap, OutMean, OutMin, OutMax,
      OutCvg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    if (csv_output != FALSE) {
      .csv.save(Alpha.Profile, OutMean, OutMin, OutMax,
        PlugIn, csv_output, CVG,
        def.fname.csv = "DP_output.csv", xx
      )
    }
    if (PlugIn == FALSE) {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(
          alpha = Alpha.Profile,
          Mean = OutMean[, l], Lower.Quantile = OutMin[, l], Upper.Quantile = OutMax[, l]
        )
      }
      names(li) <- colnames(xx)
    }
    else {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(Alpha = Alpha.Profile, PlugIn = OutMean[, l])
      }
      names(li) <- colnames(xx)
    }
    if (length(Alpha.Profile) <= 1) {
      warning("For graphics alpha must be a vector of lenght greater than 1")
    }
    y.label <- "Effective Number of Species"
    title <- "Diversity Profile:"
    single.fname <- "ENS_"
    max.y <- (max(OutMean) + 3)
    min.y <- min(OutMean) - 3

    if (pooled_graph != FALSE & length(alpha) > 1) {
      .pooled.graph(Alpha.Profile, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, pooled_graph, PlugIn,
        def.fname.pdf = "divo_ENS_AllPopulations.pdf", xx
      )
    }
    if (single_graph != FALSE & length(alpha) > 1) {
      .single.graph(
        Alpha.Profile, xx, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, single_graph, PlugIn, single.fname
      )
    }
  } else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
  return(li)
}
dp.ht <- function(x, alpha = seq(0.1, 2, 0.1), CI = 0.95, resample = 100,
                  single_graph = FALSE, pooled_graph = FALSE, csv_output = FALSE,
                  PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (0 < CI & CI < 1) {
    f <- "RE"
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    Alpha.Profile <- alpha
    if (min(Alpha.Profile) <= 0) {
      stop("alpha must be larger than 0")
    }
    DP_nAP <- length(Alpha.Profile)
    dimX <- dim(xx)
    OutMean <- mat.or.vec(DP_nAP, dimX[2])
    OutMin <- mat.or.vec(DP_nAP, dimX[2])
    OutMax <- mat.or.vec(DP_nAP, dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "DP", xx, "DP", 1, resample, CI, "HT", PlugIn, size,
      Alpha.Profile, CVG, saveBootstrap, OutMean, OutMin, OutMax,
      OutCvg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    if (csv_output != FALSE) {
      .csv.save(Alpha.Profile, OutMean, OutMin, OutMax,
        PlugIn, csv_output, CVG,
        def.fname.csv = "DP_output.csv", xx
      )
    }
    if (PlugIn == FALSE) {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(
          alpha = Alpha.Profile,
          Mean = OutMean[, l], Lower.Quantile = OutMin[, l],
          Upper.Quantile = OutMax[, l]
        )
      }
      names(li) <- colnames(xx)
    }
    else {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(Alpha = Alpha.Profile, PlugIn = OutMean[, l])
      }
      names(li) <- colnames(xx)
    }
    if (length(Alpha.Profile) <= 1) {
      warning("For graphics alpha must be a vector of lenght greater than 1")
    }
    y.label <- "Diversity Index"
    title <- "Diversity Profile with Horvitz-Thompson Adjustment:"
    single.fname <- "DP.HT_"
    max.y <- (max(OutMean) + .3)
    min.y <- min(OutMean) - .3
    if (pooled_graph != FALSE & length(alpha) > 1) {
      .pooled.graph(Alpha.Profile, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, pooled_graph, PlugIn,
        def.fname.pdf = "divo_DP.HT_AllPopulations.pdf", xx
      )
    }
    if (single_graph != FALSE & length(alpha) > 1) {
      .single.graph(
        Alpha.Profile, xx, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, single_graph, PlugIn, single.fname
      )
    }
  } else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
  return(li)
}

ens.ht <- function(x, alpha = seq(0.1, 2, 0.1), CI = 0.95, resample = 100,
                   single_graph = FALSE, pooled_graph = FALSE, csv_output = FALSE,
                   PlugIn = FALSE, size = 1, CVG = FALSE, saveBootstrap = FALSE) {
  xx <- x
  resample <- round(resample, 0)
  xx[is.na(xx)] <- 0
  if (resample < 1) {
    stop("resample must be greater than 1")
  }
  if (size < 0.1) {
    stop("size must be greater than or equal to 0.1")
  }
  if (0 < CI & CI < 1) {
    if (class(xx)[1] != "matrix") {
      xx <- as.matrix(xx)
    }
    Alpha.Profile <- alpha
    if (min(Alpha.Profile) <= 0) {
      stop("alpha must be larger than 0")
    }
    DP_nAP <- length(Alpha.Profile)
    dimX <- dim(xx)
    OutMean <- mat.or.vec(DP_nAP, dimX[2])
    OutMin <- mat.or.vec(DP_nAP, dimX[2])
    OutMax <- mat.or.vec(DP_nAP, dimX[2])
    OutCvg <- mat.or.vec(dimX[2], 1)
    .Call(
      "read", "DP", xx, "ENS", 1, resample, CI, "HT", PlugIn, size,
      Alpha.Profile, CVG, saveBootstrap, OutMean, OutMin, OutMax,
      OutCvg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    if (csv_output != FALSE) {
      .csv.save(Alpha.Profile, OutMean, OutMin, OutMax,
        PlugIn, csv_output, CVG,
        def.fname.csv = "DP_output.csv", xx
      )
    }
    if (PlugIn == FALSE) {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(
          alpha = Alpha.Profile,
          Mean = OutMean[, l], Lower.Quantile = OutMin[, l], Upper.Quantile = OutMax[, l]
        )
      }
      names(li) <- colnames(xx)
    }
    else {
      li <- list()
      for (l in 1:dim(OutMean)[2]) {
        li[[l]] <- cbind(Alpha = Alpha.Profile, PlugIn = OutMean[, l])
      }
      names(li) <- colnames(xx)
    }
    if (length(Alpha.Profile) <= 1) {
      warning("For graphics alpha must be a vector of lenght greater than 1")
    }
    y.label <- "Effective Number of Species"
    title <- "Diversity Profile with Horvitz-Thompson Adjustment:"
    single.fname <- "ENS.HT_"
    max.y <- (max(OutMean) + 3)
    min.y <- min(OutMean) - 3
    if (pooled_graph != FALSE & length(alpha) > 1) {
      .pooled.graph(Alpha.Profile, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, pooled_graph, PlugIn,
        def.fname.pdf = "divo_ENS.HT_AllPopulations.pdf", xx
      )
    }
    if (single_graph != FALSE & length(alpha) > 1) {
      .single.graph(
        Alpha.Profile, xx, OutMean, OutMax, OutMin, y.label, min.y,
        max.y, title, single_graph, PlugIn, single.fname
      )
    }
  } else {
    stop("Confidence interval must be between 0 and 1; default CI=0.95")
  }
  return(li)
}



.csv.save <- function(Alpha.Profile, OutMean, OutMin, OutMax, PlugIn, csv_output,
                      CVG, def.fname.csv, xx) {
  out <- as.data.frame(Alpha.Profile)
  if (PlugIn == TRUE) {
    out <- cbind(out, OutMean)
    colnames(out) <- c("alpha", colnames(xx))
  }
  else {
    for (l in 1:dim(OutMean)[2]) {
      tmp <- cbind(OutMean[, l], OutMin[, l], OutMax[, l])
      colnames(tmp) <- paste(as.character(colnames(xx)[l]),
        c("Mean", "Lower.Quantile", "Upper.Quantile"),
        sep = "."
      )
      out <- cbind(out, tmp)
    }
  }
  fnm <- paste(csv_output, ".csv", sep = "")
  if (csv_output == TRUE) {
    fnm <- def.fname.csv
  }
  write.csv(out, file = fnm, row.names = FALSE)
}

.pooled.graph <- function(Alpha.Profile, OutMean, OutMax, OutMin, y.label, min.y,
                          max.y, title, pooled_graph, PlugIn, def.fname.pdf, xx) {
  lcol <- rainbow(dim(OutMean)[2])
  if (PlugIn == TRUE) {
    plot(Alpha.Profile, OutMean[, 1],
      col = "blue", pch = "", ylim = c(min.y, max.y),
      main = paste(title, "All Populations"), xlab = "Order (alpha)",
      ylab = y.label, cex.axis = 0.8, cex.lab = 0.9, cex.main = 0.9
    )
    for (i in 1:dim(OutMean)[2]) {
      lines(Alpha.Profile, OutMean[, i], type = "l", col = lcol[i])
    }
  }
  else {
    plot(Alpha.Profile, OutMean[, 1],
      col = "blue", pch = "",
      ylim = c(min.y, max.y), main = paste(title, "All Populations"),
      xlab = "Order (alpha)", ylab = y.label, cex.axis = 0.8, cex.lab = 0.9,
      cex.main = 0.9
    )
    for (i in 1:dim(OutMean)[2]) {
      lines(Alpha.Profile, OutMean[, i], type = "l", col = lcol[i])
      lines(Alpha.Profile, OutMean[, i] + OutMax[, i], lty = 3, col = lcol[i])
      lines(Alpha.Profile, OutMean[, i] - OutMin[, i], lty = 3, col = lcol[i])
    }
  }
  try(legend("topright", colnames(xx),
    cex = 0.5,
    col = c(lcol[1:length(colnames(xx))]), lty = 1:1,
    lwd = 1, bty = "n", xpd = NA
  ), silent = TRUE)
  if (pooled_graph == TRUE) {
    dev.copy(pdf, def.fname.pdf)
  } else {
    dev.copy(pdf, paste(pooled_graph, ".pdf", sep = ""))
  }
  dev.off()
  return()
}


.single.graph <- function(Alpha.Profile, xx, OutMean, OutMax, OutMin, y.label,
                          min.y, max.y, title, single_graph, PlugIn,
                          single.fname) {
  if (PlugIn == TRUE) {
    for (i in 1:dim(OutMean)[2]) {
      plot(Alpha.Profile, OutMean[, i],
        type = "l", ylim = c(min.y, max.y),
        main = colnames(xx)[i], xlab = "Order (alpha)", ylab = y.label,
        cex.axis = 0.8, cex.lab = 0.9, cex.main = 0.9
      )
      if (single_graph == TRUE) {
        dev.copy(pdf, paste(single.fname,
          as.character(colnames(xx)[i]), ".pdf",
          sep = ""
        ))
      } else {
        dev.copy(pdf, paste(single_graph, "_",
          as.character(colnames(xx)[i]), ".pdf",
          sep = ""
        ))
      }
      dev.off()
    }
  } else {
    for (i in 1:dim(OutMean)[2]) {
      plot(Alpha.Profile, OutMean[, i],
        type = "l", ylim = c(min.y, max.y),
        main = colnames(xx)[i], xlab = "Order (alpha)", ylab = y.label,
        cex.axis = 0.8, cex.lab = 0.9, cex.main = 0.9
      )
      lines(Alpha.Profile, OutMean[, i] + OutMax[, i], lty = 3)
      lines(Alpha.Profile, OutMean[, i] - OutMin[, i], lty = 3)
      if (single_graph == TRUE) {
        dev.copy(pdf, paste(single.fname,
          as.character(colnames(xx)[i]), ".pdf",
          sep = ""
        ))
      } else {
        dev.copy(pdf, paste(single_graph, "_", as.character(colnames(xx)[i]),
          ".pdf",
          sep = ""
        ))
      }
      dev.off()
    }
  }
  return()
}
