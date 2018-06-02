# set as working directory
setwd("/home/danara/Documents/hierfstat/ROD");
#### Load libraries ####
library(adegenet)
library(pegas)
library(foreach)
library(doParallel)
library(spatstat)
library(tibble)
library(windowscanr)
library(dplyr)

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
          if (sum(ag %in% allel[k]) == 1) 
            h[[j]][i, k] <- h[[j]][i, k] + tmp$genotype[l]
        }
      }
    }
  }
  obj <- matrix(0, nloci, 1)
  for (i in 1:length(p)) {
    s1 <- sum(unlist(p[[i]][1,1:2]))
    a1 <- (p[[i]][1,1]/s1)^2
    b1 <- (p[[i]][1,2]/s1)^2
    di1 <- 1-(a1+b1)
    s2 <- sum(unlist(p[[i]][2,1:2]))
    a2 <- (p[[i]][2,1]/s2)^2
    b2 <- (p[[i]][2,2]/s2)^2
    di2 <- 1-(a2+b2)
    obj[i, 1] <- (di1 + 0.1)/(di2 + 0.1)
  }
  dimnames(obj) <- list(names(x)[attr(x, "locicol")], "ROD")
  return(obj)
}

# Import the table with chr and pos for each SNP, required for sliding window
POS <- read.csv("/home/danara/Documents/hierRODat/ROD/snp_chr_pos.csv", header= TRUE, stringsAsFactors = FALSE, sep = ",")
##### Import and transform the data ####
# Import the table with chr and pos for each SNP, required for sliding window
POS <- read.csv("/home/danara/Documents/hierfstat/Fst/snp_chr_pos.csv", header= TRUE, stringsAsFactors = FALSE, sep = ",")

# DWL
DWL <- read.table("/home/danara/Documents/hierfstat/Fst-4-pops/ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DWL_ROD.txt", sep = "\t", dec = ".", h = T,comment.char = "?")
for(i in 12:ncol(DWL)){
  DWL[, i] <-  gsub("^AA$", "A-A", DWL[, i])
  DWL[, i] <-  gsub("^TT$", "B-B", DWL[, i])
}
DWLt <- as.data.frame(t(DWL[, 12:ncol(DWL)]))
names(DWLt) <- DWL[, 1]
DWLt <- cbind(rep("DWL", nrow(DWLt)),DWLt)
names(DWLt)[1] <- "population"

# DWC
DWC <- read.table("/home/danara/Documents/hierfstat/Fst-4-pops/ld099/180422_SvevoDiversity_SNPfiltered_file_all_17340K_SNP_1765_genot_NOT_IMPUTED_LD099_maxNN_025 V2_DWC_ROD.txt", sep = "\t", dec = ".", h = T,comment.char = "?")
for(i in 12:ncol(DWC)){
  DWC[, i] <-  gsub("^AA$", "A-A", DWC[, i])
  DWC[, i] <-  gsub("^TT$", "B-B", DWC[, i])
}
DWCt <- as.data.frame(t(DWC[, 12:ncol(DWC)]))
names(DWCt) <- DWC[, 1]
DWCt <- cbind(rep("DWC", nrow(DWCt)),DWCt)
names(DWCt)[1] <- "population"

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
  ROD_unibo(tmp, pop = 1)
}
write.table(results, file="ROD_1M_permutations_DWL-DWC.txt", sep=",", row.names=T, col.names = NA, quote = FALSE)

# Prepare an empty matrix for quantiles 
N <- matrix(NA, ncol=2, nrow=nrow(results))
# calculate the distribution, quantile
snps <- ncol(loci)-1
for(j in 1:nrow(results)){
  d <- density(results[j,])
  q <- quantile(d,probs=c((1-(0.05/snps)),(1-(0.05/9946))), na.rm=TRUE)
  N[j,]=q
}

# calculate ROD effettivo senza sliding window
RODcomput <- ROD_unibo(loci)
final_N <- cbind(N, RODcomput)
final_N <- as.data.frame(final_N)
final_N$delta_ROD1 <- (final_N$ROD - final_N$V1)
final_N$delta_ROD2 <- (final_N$ROD - final_N$V2)
names(final_N) <- c("soglia.1", "soglia.2", "ROD", "delta_ROD1", "delta_ROD2")
write.table(final_N, file="ROD_permutation_threshold_DWL-DWC.txt", sep="\t", col.names=NA, row.names=T, quote = FALSE)

# calculate ROD effettivo con sliding window da 2 Mb e step 1 Mb
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

# calculate ROD effettivo
#RODcomput <- ROD_unibo(loci)
RODcomput <- as.data.frame(RODcomput)
RODcomput <- add_column(RODcomput, group = NA, pos =  NA, .before = "ROD")
RODcomput$group <- POS$chrom[match(row.names(RODcomput), POS$SNP_ID)]
RODcomput$pos <- POS$pos[match(row.names(RODcomput), POS$SNP_ID)]
ROD_pos_win <- winScan(x = RODcomput,
                       groups = "group",
                       position = "pos",
                       values = "ROD",
                       win_size = 2000000,
                       win_step = 1000000,
                       funs = "mean",
                       cores = 7)
row.has.na <- apply(ROD_pos_win, 1, function(x){any(is.na(x))})
ROD_pos_win <- ROD_pos_win[!row.has.na,]
row.has.na <- apply(pos_win, 1, function(x){any(is.na(x))})
pos_win <- pos_win[!row.has.na,]
final_ROD <- cbind(pos_win, ROD_pos_win$ROD_mean)
final_N <- as.data.frame(final_N)
final_ROD$delta_ROD1 <- (final_ROD$`ROD_pos_win$ROD_mean` - final_ROD$soglia.1)
final_ROD$delta_ROD2 <- (final_ROD$`ROD_pos_win$ROD_mean` - final_ROD$soglia.2)
names(final_ROD)[7] <- "ROD"
write.table(final_ROD, file="ROD_permutation_threshold_SW_DWL-DWC.txt", sep="\t", row.names=T, col.names = NA, quote = FALSE)
stopCluster(cl)