
# load common functions
source("common.r")

#### (1) ROD for cross-population DEW - DWL ####
# setup parallel backend to use many processors
# cores are defined in commons.r
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# number of permutations are in common.r

# load data
DEWt <- read_DEWt()
DWLt <- read_DWLt()

# bind populations by row
DD <- rbind(DEWt, DWLt)

# remove unnecessary data
rm(DWLt)
rm(DEWt)

# calc genind
all_genind <- df2genind(DD[2:ncol(DD)], NA.char = "NN", ploidy = 2, pop = DD$population, sep = '-')

# remove unnecessary data
rm(DD)

# transform loci
loci <- as.loci(all_genind)

# Get rid of monomorphic SNPs
loci <- loci[, !apply(loci, 2, function(x) length(levels(as.factor(x))) == 1)]

# SUBSET TO TEST IF IT WORKS
# loci <- loci[,1: 20,]

# setting seed
set.seed(100)

tic("Permutations")
results <- foreach(i=1:permutations, .combine=cbind, .packages = "pegas")  %dopar% {
  tmp <- loci
  tmp$population <- tmp$population[sample(1:length(tmp$population))]
  ROD_unibo(tmp, pop = 1)
}
toc()
# write.table(results, file=paste("ROD", permutations, "permutations_DEW-DWL.txt", sep="_"), sep=",", row.names=T, col.names = NA, quote = FALSE)
dir.create(file.path(current_dir, "RData"), showWarnings = FALSE)
saveRDS(results, file=paste("RData/ROD", permutations, "permutations_DEW-DWL.rds", sep="_"))

# Prepare an empty matrix for quantiles
N <- matrix(NA, ncol=3, nrow=nrow(results))

# calculate the distribution, quantile 97.5 % and 99 %
snps <- ncol(loci)-1
for(j in 1:nrow(results)){
  d <- density(results[j,])
  q <- quantile(d,probs=c((1-(0.05/snps)),(1-(0.05/9946)),(1-(0.05/5775))), na.rm=TRUE)
  N[j,]=q
}

# calculate ROD effettivo senza sliding window
RODcomput <- ROD_unibo(loci)
final_N <- cbind(N, RODcomput)
final_N <- as.data.frame(final_N)
final_N$delta_ROD1 <- (final_N$ROD - final_N$V1)
final_N$delta_ROD2 <- (final_N$ROD - final_N$V2)
final_N$delta_ROD3 <- (final_N$ROD - final_N$V3)
names(final_N) <- c("soglia.1", "soglia.2","soglia.3", "ROD", "delta_ROD1", "delta_ROD2", "delta_ROD3")
write.table(final_N, file="ROD_permutation_threshold_DEW-DWL.txt", sep="\t", col.names=NA, row.names=T, quote = FALSE)

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
names(pos_win)[5:7] <- c("soglia.1", "soglia.2", "soglia.3")

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
                       cores = cores[1]-1)
row.has.na <- apply(ROD_pos_win, 1, function(x){any(is.na(x))})
ROD_pos_win <- ROD_pos_win[!row.has.na,]
row.has.na <- apply(pos_win, 1, function(x){any(is.na(x))})
pos_win <- pos_win[!row.has.na,]
final_ROD <- cbind(pos_win, ROD_pos_win$ROD_mean)
final_N <- as.data.frame(final_N)
final_ROD$delta_ROD1 <- (final_ROD$`ROD_pos_win$ROD_mean` - final_ROD$soglia.1)
final_ROD$delta_ROD2 <- (final_ROD$`ROD_pos_win$ROD_mean` - final_ROD$soglia.2)
final_ROD$delta_ROD3 <- (final_ROD$`ROD_pos_win$ROD_mean` - final_ROD$soglia.3)
names(final_ROD)[8] <- "ROD"
write.table(final_ROD, file="ROD_permutation_threshold_SW_DEW-DWL.txt", sep="\t", row.names=T, col.names = NA, quote = FALSE)
stopCluster(cl)
