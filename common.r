
# set as working directory
# setwd("/home/danara/Documents/hierfstat/");

# get current directory as basedir
current_dir <- getwd()

# get current number of CPUs
cores <- strtoi(Sys.getenv("CORES",unset = 4))

#### Load libraries ####
library(adegenet)
library(pegas)
library(foreach)
library(doParallel)
library(spatstat)
library(tibble)
library(windowscanr)
library(dplyr)
library(tictoc)

#### Fst function from pegas package modified  ####
Fst_unibo <- function (x, pop = NULL)
{
  if (any(getPloidy(x) != 2))
    stop("Fst() requires diploid data")

  if (is.null(pop)) {
    pop <- x$population
    if (is.null(pop))
      stop("no 'population' column in x")
  }
  else {
    pop <- if (is.numeric(pop))
      x[, pop]
    else factor(pop)
  }

  r <- length(attr(pop, "levels"))
  pop <- as.integer(pop)
  nloci <- length(attr(x, "locicol"))
  ALLELES <- getAlleles(x)
  p <- vector("list", nloci)

  for (j in 1:nloci) p[[j]] <- matrix(0, r, length(ALLELES[[j]]))
  h <- p
  for (i in 1:r) {
    s <- summary(x[pop == i, ])
    for (j in 1:nloci) {
      tmp <- s[[j]]
      p[[j]][i, ] <- tmp$allele
      allel <- names(tmp$allele)
      genot <- names(tmp$genotype)
      for (k in seq_along(allel)) {
        for (l in seq_along(genot)) {
          ag <- unlist(strsplit(genot[l], "/"))
          if (sum(ag %in% allel[k]) == 1)
            h[[j]][i, k] <- h[[j]][i, k] + tmp$genotype[l]
        }
      }
    }
  }
  obj <- matrix(0, nloci, 1)
  for (j in 1:nloci) {
    nBYpop <- rowSums(p[[j]])/2
    N <- (sum(nBYpop))
    nbar <- N/r
    nC <- (N - sum(nBYpop^2)/N)/(r - 1)
    ptild <- p[[j]]/(2 * nBYpop)
    pbar <- colSums(p[[j]])/(2 * N)
    s2 <- colSums(nBYpop * (ptild - rep(pbar, each = r))^2)/((r -
                                                                1) * nbar)
    hbar <- colSums(h[[j]])/N
    A <- pbar * (1 - pbar) - (r - 1) * s2/r
    a <- nbar * (s2 - (A - hbar/4)/(nbar - 1))/nC
    b <- nbar * (A - (2 * nbar - 1) * hbar/(4 * nbar))/(nbar -
                                                          1)
    c <- hbar/2
    obj[j, 1] <- sum(a)/sum(a + b + c)
  }
  dimnames(obj) <- list(names(x)[attr(x, "locicol")], "Fst")
  obj
}



#### Load the ROD function ####
ROD_unibo <- function (x, pop = NULL) {
  pop <- as.factor(x$population)
  r <- length(attr(pop, "levels"))
  pop <- as.integer(pop)
  nloci <- length(attr(x, "locicol"))
  ALLELES <- getAlleles(x)
  p <- vector("list", nloci)
  for (j in 1:nloci) p[[j]] <- matrix(0, r, length(ALLELES[[j]]))
  h <- p
  for (i in 1:r) { #i = populations
    s <- summary(x[pop == i, ])
    for (j in 1:nloci) { # for wach marker
      tmp <- s[[j]]
      p[[j]][i, ] <- tmp$allele
      allel <- names(tmp$allele)
      genot <- names(tmp$genotype)
      for (k in seq_along(allel)) {
        for (l in (seq_along(genot))) {
          ag <- unlist(strsplit(genot[l], "/"))
          #if (sum(ag %in% allel[k]) == 1)
          #  h[[j]][i, k] <- h[[j]][i, k] + tmp$genotype[l]
        }
      }
    }
  }
  obj <- matrix(0, nloci, 1)
  for (j in 1:nloci) {
    nBYpop <- rowSums(p[[j]])
    a_b <- (p[[j]]/nBYpop)^2
    di <- 1-rowSums(a_b)
    obj[j, 1] <- (di[1] + 0.1)/(di[2] + 0.1)
  }
  dimnames(obj) <- list(names(x)[attr(x, "locicol")], "ROD")
  return(obj)
}


##### Import and transform the data ####
# Import the table with chr and pos for each SNP, required for sliding window
POS <- read.csv(file.path(current_dir, "ld099/snp_chr_pos.csv"), header= TRUE, stringsAsFactors = FALSE, sep = ",", nrow=1000)

read_DEWt <- function () {
  # DEW
  DEW <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DEW_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?", nrows=1000)
  for(i in 12:ncol(DEW)){
    DEW[, i] <-  gsub("^AA$", "A-A", DEW[, i])
    DEW[, i] <-  gsub("^TT$", "B-B", DEW[, i])
  }
  DEWt <- as.data.frame(t(DEW[, 12:ncol(DEW)]))
  names(DEWt) <- DEW[, 1]
  DEWt <- cbind(rep("DEW", nrow(DEWt)),DEWt)
  names(DEWt)[1] <- "population"

  return(DEWt)
}

read_DWLt <- function () {
  # DWL
  DWL <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DWL_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?", nrows=1000)
  for(i in 12:ncol(DWL)){
    DWL[, i] <-  gsub("^AA$", "A-A", DWL[, i])
    DWL[, i] <-  gsub("^TT$", "B-B", DWL[, i])
  }
  DWLt <- as.data.frame(t(DWL[, 12:ncol(DWL)]))
  names(DWLt) <- DWL[, 1]
  DWLt <- cbind(rep("DWL", nrow(DWLt)),DWLt)
  names(DWLt)[1] <- "population"

  return(DWLt)
}

read_DWCt <- function() {
  # DWC
  DWC <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DWC_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?", nrows=1000)
  for(i in 12:ncol(DWC)){
    DWC[, i] <-  gsub("^AA$", "A-A", DWC[, i])
    DWC[, i] <-  gsub("^TT$", "B-B", DWC[, i])
  }
  DWCt <- as.data.frame(t(DWC[, 12:ncol(DWC)]))
  names(DWCt) <- DWC[, 1]
  DWCt <- cbind(rep("DWC", nrow(DWCt)),DWCt)
  names(DWCt)[1] <- "population"

  return(DWCt)
}

read_WEWt <- function() {
  # WEW
  WEW <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_WEW_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?", nrows = 1000)
  for(i in 12:ncol(WEW)){
    WEW[, i] <-  gsub("^AA$", "A-A", WEW[, i])
    WEW[, i] <-  gsub("^TT$", "B-B", WEW[, i])
  }
  WEWt <- as.data.frame(t(WEW[, 12:ncol(WEW)]))
  names(WEWt) <- WEW[, 1]
  WEWt <- cbind(rep("WEW", nrow(WEWt)),WEWt)
  names(WEWt)[1] <- "population"

  return(WEWt)
}
