# set as working directory
# setwd("/home/danara/Documents/hierfstat/");

# get current directory as basedir
current_dir <- getwd()

#### Load libraries ####
library(adegenet)
library(pegas)
library(foreach)
library(doParallel)
library(spatstat)
library(tibble)
library(windowscanr)
library(dplyr)

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

##### Import and transform the data ####
# Import the table with chr and pos for each SNP, required for sliding window
POS <- read.csv(file.path(current_dir, "ld099/snp_chr_pos.csv"), header= TRUE, stringsAsFactors = FALSE, sep = ",")
# DEW
DEW <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DEW_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?")
for(i in 12:ncol(DEW)){
  DEW[, i] <-  gsub("^AA$", "A-A", DEW[, i])
  DEW[, i] <-  gsub("^TT$", "B-B", DEW[, i])
}
DEWt <- as.data.frame(t(DEW[, 12:ncol(DEW)]))
names(DEWt) <- DEW[, 1]
DEWt <- cbind(rep("DEW", nrow(DEWt)),DEWt)
names(DEWt)[1] <- "population"

# DWL
DWL <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DWL_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?")
for(i in 12:ncol(DWL)){
  DWL[, i] <-  gsub("^AA$", "A-A", DWL[, i])
  DWL[, i] <-  gsub("^TT$", "B-B", DWL[, i])
}
DWLt <- as.data.frame(t(DWL[, 12:ncol(DWL)]))
names(DWLt) <- DWL[, 1]
DWLt <- cbind(rep("DWL", nrow(DWLt)),DWLt)
names(DWLt)[1] <- "population"

# DWC
DWC <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DWC_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?")
for(i in 12:ncol(DWC)){
  DWC[, i] <-  gsub("^AA$", "A-A", DWC[, i])
  DWC[, i] <-  gsub("^TT$", "B-B", DWC[, i])
}
DWCt <- as.data.frame(t(DWC[, 12:ncol(DWC)]))
names(DWCt) <- DWC[, 1]
DWCt <- cbind(rep("DWC", nrow(DWCt)),DWCt)
names(DWCt)[1] <- "population"

# WEW
WEW <- read.table(file.path(current_dir, "ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_WEW_ROD.txt"), sep = "\t", dec = ".", h = T,comment.char = "?")
for(i in 12:ncol(WEW)){
  WEW[, i] <-  gsub("^AA$", "A-A", WEW[, i])
  WEW[, i] <-  gsub("^TT$", "B-B", WEW[, i])
}
WEWt <- as.data.frame(t(WEW[, 12:ncol(WEW)]))
names(WEWt) <- WEW[, 1]
WEWt <- cbind(rep("WEW", nrow(WEWt)),WEWt)
names(WEWt)[1] <- "population"

#### (1) Fst for cross-population DEW - DWL ####
# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
#registerDoParallel(cores=2)
registerDoParallel(cl)
DD <- rbind(DEWt, DWLt)
all_genind <- df2genind(DD[2:ncol(DD)], NA.char = "NN", ploidy = 2, pop = DD$population, sep = "-")
loci <- as.loci(all_genind)
loci <- loci[, !apply(loci, 2, function(x) length(levels(as.factor(x))) == 1)] # Get rid of monomorphic SNPs
#loci <- loci[,1: 20,] #SUBSET TO TEST IF IT WORKS
set.seed(100)
results <- foreach(i=1:3, .combine=cbind, .packages = "pegas")  %dopar% {
  tmp <- loci
  tmp$population <- tmp$population[sample(1:length(tmp$population))]
  Fst_unibo(tmp, pop = 1)
}
write.table(results, file="Fst_1M_permutations_DEW-DWL.txt", sep=",", row.names=T, col.names = NA, quote = FALSE)

# Prepare an empty matrix for quantiles
N <- matrix(NA, ncol=2, nrow=nrow(results))
# calculate the distribution, quantile 97.5 % and 99 %
snps <- ncol(loci)-1
for(j in 1:nrow(results)){
  d <- density(results[j,])
  q <- quantile(d,probs=c((1-(0.05/snps)),(1-(0.05/9946))), na.rm=TRUE)
  N[j,]=q
}

# calculate Fst effettivo senza sliding window
Fstcomput <- Fst_unibo(loci)
final_N <- cbind(N, Fstcomput)
final_N <- as.data.frame(final_N)
final_N$delta_Fst1 <- (final_N$Fst - final_N$V1)
final_N$delta_Fst2 <- (final_N$Fst - final_N$V2)
names(final_N) <- c("soglia.1", "soglia.2", "Fst", "delta_Fst1", "delta_Fst2")
write.table(final_N, file="Fst_permutation_threshold_DEW-DWL.txt", sep="\t", col.names=NA, row.names=T, quote = FALSE)

# calculate Fst effettivo con sliding window da 2 Mb e step 1 Mb
row.names(N) <- row.names(results)
results2 <- as.data.frame(N)
# Position sliding window for quantiles
results2 <- add_column(results2, group = NA, pos =  NA, .before = "V1")
results2$group <- POS$chrom[match(row.names(results2), POS$SNP_ID)]
results2$pos <- POS$pos[match(row.names(results2), POS$SNP_ID)]

pos_win <- foreach(i = 3:ncol(results2), .combine=cbind, .packages = c("windowscanr", "dplyr")) %dopar% {
  pos_win_tmp <- winScan(x = results2,
                         groups = "group",
                         position = "pos",
                         values = as.character(names(results2)[i]),
                         win_size = 2000000,
                         win_step = 1000000,
                         funs = "mean")
  if(i == 3){pos_win_tmp[, c(1:4, 6)]}else{pos_win_tmp[, 6]}
}
names(pos_win)[5:6] <- c("soglia.1", "soglia.2")

# calculate Fst effettivo
#Fstcomput <- Fst_unibo(loci)
Fstcomput <- as.data.frame(Fstcomput)
Fstcomput <- add_column(Fstcomput, group = NA, pos =  NA, .before = "Fst")
Fstcomput$group <- POS$chrom[match(row.names(Fstcomput), POS$SNP_ID)]
Fstcomput$pos <- POS$pos[match(row.names(Fstcomput), POS$SNP_ID)]
Fst_pos_win <- winScan(x = Fstcomput,
                       groups = "group",
                       position = "pos",
                       values = "Fst",
                       #values = c(results2[, 3: ncol(results2)]),
                       win_size = 2000000,
                       win_step = 1000000,
                       funs = "mean",
                       cores = 7)
row.has.na <- apply(Fst_pos_win, 1, function(x){any(is.na(x))})
Fst_pos_win <- Fst_pos_win[!row.has.na,]
row.has.na <- apply(pos_win, 1, function(x){any(is.na(x))})
pos_win <- pos_win[!row.has.na,]
final_Fst <- cbind(pos_win, Fst_pos_win$Fst_mean)
final_N <- as.data.frame(final_N)
final_Fst$delta_Fst1 <- (final_Fst$`Fst_pos_win$Fst_mean` - final_Fst$soglia.1)
final_Fst$delta_Fst2 <- (final_Fst$`Fst_pos_win$Fst_mean` - final_Fst$soglia.2)
names(final_Fst)[7] <- "Fst"
write.table(final_Fst, file="Fst_permutation_threshold_SW_DEW-DWL.txt", sep="\t", row.names=T, col.names = NA, quote = FALSE)
stopCluster(cl)

#### (2) Fst for cross-population WEW - DEW ####
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
DD <- rbind(WEWt, DEWt)
all_genind <- df2genind(DD[2:ncol(DD)], NA.char = "NN", ploidy = 2, pop = DD$population, sep = "-")
loci <- as.loci(all_genind)
loci <- loci[, !apply(loci, 2, function(x) length(levels(as.factor(x))) == 1)] # Get rid of monomorphic SNPs
#loci <- loci[,1: 20,] #SUBSET TO TEST IF IT WORKS
set.seed(100)
results <- foreach(i=1:1000000, .combine=cbind, .packages = "pegas")  %dopar% {
  tmp <- loci
  tmp$population <- tmp$population[sample(1:length(tmp$population))]
  Fst_unibo(tmp, pop = 1)
}
write.table(results, file="Fst_1M_permutations_WEW-DEW.txt", sep=",", row.names=T, col.names = NA, quote = FALSE)

# Prepare an empty matrix for quantiles
N <- matrix(NA, ncol=2, nrow=nrow(results))
# calculate the distribution, quantile
snps <- ncol(loci)-1
for(j in 1:nrow(results)){
  d <- density(results[j,])
  q <- quantile(d,probs=c((1-(0.05/snps)),(1-(0.05/9946))), na.rm=TRUE)
  N[j,]=q
}

# calculate Fst effettivo senza sliding window
Fstcomput <- Fst_unibo(loci)
final_N <- cbind(N, Fstcomput)
final_N <- as.data.frame(final_N)
final_N$delta_Fst1 <- (final_N$Fst - final_N$V1)
final_N$delta_Fst2 <- (final_N$Fst - final_N$V2)
names(final_N) <- c("soglia.1", "soglia.2", "Fst", "delta_Fst1", "delta_Fst2")
write.table(final_N, file="Fst_permutation_threshold_WEW-DEW.txt", sep="\t", col.names=NA, row.names=T, quote = FALSE)

# calculate Fst effettivo con sliding window da 2 Mb e step 1 Mb
row.names(N) <- row.names(results)
results2 <- as.data.frame(N)
# Position sliding window for quantiles
results2 <- add_column(results2, group = NA, pos =  NA, .before = "V1")
results2$group <- POS$chrom[match(row.names(results2), POS$SNP_ID)]
results2$pos <- POS$pos[match(row.names(results2), POS$SNP_ID)]

pos_win <- foreach(i = 3:ncol(results2), .combine=cbind, .packages = c("windowscanr", "dplyr")) %dopar% {
  pos_win_tmp <- winScan(x = results2,
                         groups = "group",
                         position = "pos",
                         values = as.character(names(results2)[i]),
                         win_size = 2000000,
                         win_step = 1000000,
                         funs = "mean")
  if(i == 3){pos_win_tmp[, c(1:4, 6)]}else{pos_win_tmp[, 6]}
}
names(pos_win)[5:6] <- c("soglia.1", "soglia.2")

# calculate Fst effettivo
#Fstcomput <- Fst_unibo(loci)
Fstcomput <- as.data.frame(Fstcomput)
Fstcomput <- add_column(Fstcomput, group = NA, pos =  NA, .before = "Fst")
Fstcomput$group <- POS$chrom[match(row.names(Fstcomput), POS$SNP_ID)]
Fstcomput$pos <- POS$pos[match(row.names(Fstcomput), POS$SNP_ID)]
Fst_pos_win <- winScan(x = Fstcomput,
                       groups = "group",
                       position = "pos",
                       values = "Fst",
                       #values = c(results2[, 3: ncol(results2)]),
                       win_size = 2000000,
                       win_step = 1000000,
                       funs = "mean",
                       cores = 7)
row.has.na <- apply(Fst_pos_win, 1, function(x){any(is.na(x))})
Fst_pos_win <- Fst_pos_win[!row.has.na,]
row.has.na <- apply(pos_win, 1, function(x){any(is.na(x))})
pos_win <- pos_win[!row.has.na,]
final_Fst <- cbind(pos_win, Fst_pos_win$Fst_mean)
final_N <- as.data.frame(final_N)
final_Fst$delta_Fst1 <- (final_Fst$`Fst_pos_win$Fst_mean` - final_Fst$soglia.1)
final_Fst$delta_Fst2 <- (final_Fst$`Fst_pos_win$Fst_mean` - final_Fst$soglia.2)
names(final_Fst)[7] <- "Fst"
write.table(final_Fst, file="Fst_permutation_threshold_SW_WEW-DEW.txt", sep="\t", row.names=T, col.names = NA, quote = FALSE)
stopCluster(cl)

#### (3) Cross-population DWL - DWC ####
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
DD <- rbind(DWLt, DWCt)
all_genind <- df2genind(DD[2:ncol(DD)], NA.char = "NN", ploidy = 2, pop = DD$population, sep = "-")
loci <- as.loci(all_genind)
loci <- loci[, !apply(loci, 2, function(x) length(levels(as.factor(x))) == 1)] # Get rid of monomorphic SNPs
#loci <- loci[,1: 20,] #SUBSET TO TEST IF IT WORKS
set.seed(100)
results <- foreach(i=1:1000000, .combine=cbind, .packages = "pegas")  %dopar% {
  tmp <- loci
  tmp$population <- tmp$population[sample(1:length(tmp$population))]
  Fst_unibo(tmp, pop = 1)
}
write.table(results, file="Fst_1M_permutations_DWL-DWC.txt", sep=",", row.names=T, col.names = NA, quote = FALSE)

# Prepare an empty matrix for quantiles
N <- matrix(NA, ncol=2, nrow=nrow(results))
# calculate the distribution, quantile
snps <- ncol(loci)-1
for(j in 1:nrow(results)){
  d <- density(results[j,])
  q <- quantile(d,probs=c((1-(0.05/snps)),(1-(0.05/9946))), na.rm=TRUE)
  N[j,]=q
}

# calculate Fst effettivo senza sliding window
Fstcomput <- Fst_unibo(loci)
final_N <- cbind(N, Fstcomput)
final_N <- as.data.frame(final_N)
final_N$delta_Fst1 <- (final_N$Fst - final_N$V1)
final_N$delta_Fst2 <- (final_N$Fst - final_N$V2)
names(final_N) <- c("soglia.1", "soglia.2", "Fst", "delta_Fst1", "delta_Fst2")
write.table(final_N, file="Fst_permutation_threshold_DWL-DWC.txt", sep="\t", col.names=NA, row.names=T, quote = FALSE)

# calculate Fst effettivo con sliding window da 2 Mb e step 1 Mb
row.names(N) <- row.names(results)
results2 <- as.data.frame(N)
# Position sliding window for quantiles
results2 <- add_column(results2, group = NA, pos =  NA, .before = "V1")
results2$group <- POS$chrom[match(row.names(results2), POS$SNP_ID)]
results2$pos <- POS$pos[match(row.names(results2), POS$SNP_ID)]

pos_win <- foreach(i = 3:ncol(results2), .combine=cbind, .packages = c("windowscanr", "dplyr")) %dopar% {
  pos_win_tmp <- winScan(x = results2,
                         groups = "group",
                         position = "pos",
                         values = as.character(names(results2)[i]),
                         win_size = 2000000,
                         win_step = 1000000,
                         funs = "mean")
  if(i == 3){pos_win_tmp[, c(1:4, 6)]}else{pos_win_tmp[, 6]}
}
names(pos_win)[5:6] <- c("soglia.1", "soglia.2")

# calculate Fst effettivo
#Fstcomput <- Fst_unibo(loci)
Fstcomput <- as.data.frame(Fstcomput)
Fstcomput <- add_column(Fstcomput, group = NA, pos =  NA, .before = "Fst")
Fstcomput$group <- POS$chrom[match(row.names(Fstcomput), POS$SNP_ID)]
Fstcomput$pos <- POS$pos[match(row.names(Fstcomput), POS$SNP_ID)]
Fst_pos_win <- winScan(x = Fstcomput,
                       groups = "group",
                       position = "pos",
                       values = "Fst",
                       win_size = 2000000,
                       win_step = 1000000,
                       funs = "mean",
                       cores = 7)
row.has.na <- apply(Fst_pos_win, 1, function(x){any(is.na(x))})
Fst_pos_win <- Fst_pos_win[!row.has.na,]
row.has.na <- apply(pos_win, 1, function(x){any(is.na(x))})
pos_win <- pos_win[!row.has.na,]
final_Fst <- cbind(pos_win, Fst_pos_win$Fst_mean)
final_N <- as.data.frame(final_N)
final_Fst$delta_Fst1 <- (final_Fst$`Fst_pos_win$Fst_mean` - final_Fst$soglia.1)
final_Fst$delta_Fst2 <- (final_Fst$`Fst_pos_win$Fst_mean` - final_Fst$soglia.2)
names(final_Fst)[7] <- "Fst"
write.table(final_Fst, file="Fst_permutation_threshold_SW_DWL-DWC.txt", sep="\t", row.names=T, col.names = NA, quote = FALSE)
stopCluster(cl)
