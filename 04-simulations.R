##%######################################################%##
#                                                          #
####               we perform the various               ####
####           simulations used in this study           ####
#                                                          #
##%######################################################%##


source("~/Mm/WZ_functions.R")

setwd("~/Documents/scrap/Mm/simu/")
setwd("~/Mm/simu")

# we import the data we need for the simulations:
# data.tables listing positions of SNPs, contig, and read pair identifiers at those positions
files <- list.files(pattern = "reads.RDS")
F1reads <- setNames(lapply(files, readRDS), splitToColumns(files, ".", 1))

# we place the heterogametic sex first
F1reads = F1reads[c(3,1,4,2)]

# we optimize speed by using vectorization as much as possible. All contigs are treated at once.
# for this , we create a table listing the read pairs for each contig
# (without duplicates) as a read pair spanning several SNPs appears multiple times
uReads <- lapply(F1reads, function(x) x[!duplicated(readPair), .(readPair, contig)])

# we also create vectors of contig identifiers corresponding to each SNP in each family
orderedContigs <- lapply(F1reads[c(1, 3)], function(x) {

  # we first obtain the last SNP position per contig
  dt <- x[, max(pos), by = contig]

  # we order contigs according to ascending SNP positions
  setorder(dt, V1)

  # we compute the number of SNPs per contig
  nSNPs <- dt[, c(V1[1], diff(V1))]

  # and return an integer vector, such that vector[x] returns the contig for SNP x
  # this help simulation speed
  rep(dt$contig, nSNPs)
})

# we compute the number of SNPs per family
nSNPs <- lengths(orderedContigs)
names(nSNPs) <- splitToColumns(names(nSNPs), "_", 2)

# and create of vector for all contigs identifiers, which will be useful
uContigs <- sort(unique(unlist(orderedContigs, use.names = F)))
nCont = max(uContigs)

ni = setNames(c(16L, 20L, 30L, 30L), names(F1reads)) # the number of individuals per pool in our experiment (in the same order as the pools)

simulate <- function(d = 0.5) {
  # this is the function that computes nrec given a simulated distance d to the
  # SDR in Morgans. It does it for all the contigs and all pools at once (it
  # cannot do it for a single contig or pool)
  
  # d can be a vector, in which case contigs can have different values of d
  d <- rep(d, length.out = nCont)

  # based on d, we sample the numbers of chromosomes carrying haplotype Z (or X) in pools
  nZ <- as.data.table(lapply(ni, function(x) rbinom(nCont, x, d)))

  # we need to correct these for the homogametic sex
  poolNames = names(nZ)
  nZ[, poolNames[c(2, 4)] := .(ni[2] - get(poolNames[2]), ni[4] - get(poolNames[4]))]
  
  # so for each pool, nZ[x, pool] gives the number of
  # chromosomes carrying haplotype Z for the contig in the pool

  # for each family we sample the haplotype carrying each maternal allele (= SNP type)
  # 1 = W, 2 = Z.
  haplo <- lapply(nSNPs, function(x) sample(1:2, x, replace = T))

  # we replicate it for each pool of the family
  haplo <- rep(haplo, each = 2)
  # so haplo[[pool]][x] returns the type of SNP x in pool pool

  # we determine whether a read carries the maternal allele at each SNP according to nZ
  for (pool in 1:length(F1reads)) {
    
    # we retrieve the number of chromosomes carrying haplotype z in the pool
    # it's easier to turn the nz data table into a list for this
    nZs <- as.list(nZ)[[pool]]
    
    # We attribute each read to maternal DNA (T) or paternal DNA (F) with same probability
    readHapl <- sample(c(T, F), nrow(uReads[[pool]]), T)
    # so readHapl[x] tells whether read x comes from maternal DNA

    # we assign maternal reads to Z- or W-linked haplotypes. To do so, we use a binomial
    # distribution with probability = the frequency of chromosomes carrying
    # haplotype Z among the ni maternally inherited chromosomes. +1L is used to
    # convert 0,1 output by rbinom() into 1,2 (1 = W, 2 = Z). paternal reads
    # will have number 0
    readHapl[readHapl] <- uReads[[pool]][readHapl, rbinom(.N, 1L, nZs[contig] / ni[pool]) + 1L]

    # the read carries the maternal allele if the source haplotype of the
    # maternal allele is the same as that of the read
    F1reads[[pool]][, maternal := readHapl[readPair] == haplo[[pool]][pos]]
  }


  # to infer SNP type, we get the number of reads for each allele, per pool
  getCounts <- function(f1reads, snps) {
    counts <- f1reads[, tabulate(pos * 2L - 1L + !maternal, nbins = snps * 2L)]
    matrix(counts, ncol = 2, byrow = T)
  }

  counts <- Map(getCounts, F1reads, rep(nSNPs, each = 2))

  # we split the results per family, as SNP type is relevant to a family
  counts <- split(counts, splitToColumns(names(counts), "_", 2))

  getType <- function(twoPools) {
    # we combine matrices of the two pools in a data.table
    twoPools <- do.call(data.table, twoPools)
    setnames(twoPools, stri_c("V", 1:4))
    type <- twoPools[, 2L - (V1 / (V1 + V2) > V3 / (V3 + V4))]
    w = which(is.na(type))
    type[w] = w %% 2L + 1L
    type
  }

  types <- lapply(counts, getType)

  # we get the per-contig proportion of SNPs for which SNP type was correctly inferred
  # which represents the accuracy of haplotype phasing
  prop <- Map(function(hapl, type, contig) {
    dt <- data.table(test = hapl == type, contig)
    dt <- dt[, .(prop = mean(test, na.rm = T)), by = contig]
    prop <- dt[match(uContigs, contig), -1, with = F]
  }, 
    haplo[c(1, 3)], types, orderedContigs)

  # we compute the f posteriors
  computeProb <- function(f1reads, group, n) {
    # we make the counts of read pairs covering the inferred haplotypes in the
    # same way we did for real data
    readCount <- f1reads[, .(contig = 1:nCont, 
                             count = tabulate(contig[!duplicated(readPair)], nbins = nCont)), 
                         by = .(group = group[pos], maternal)]
    readCount <- dcast(readCount, contig ~ ifelse(maternal, "c", "r") + group, 
                       value.var = "count", sep = "", fill = 0L)

    # if some cases, the c2 column may be missing (if nZ = 0 for all contigs)
    if (ncol(readCount) < 5L) readCount[, c2 := 0L]
    
    setcolorder(readCount, c("contig", "c1", "r1", "c2", "r2"))
    readCount[, c("r1", "r2") := .(c1 + r1, c2 + r2)]
    names(n) <- NULL
    probas <- do.call(pFgivenCountsHapl, c(readCount[uContigs, -"contig"], ni=n))
  }

  # we compute probabilities with inferred SNP type (phased haplotypes)
  probasInferred <- Map(computeProb, 
                        F1reads, 
                        list(types[[1]], 3L - types[[1]], types[[2]], 3L - types[[2]]), ni)

  # and with simulated SNP type (true haplotype)
  probasSimulated <- Map(computeProb, 
                         F1reads, 
                         list(haplo[[1]], 3L - haplo[[2]], haplo[[3]], 3L - haplo[[4]]), ni)

  computeNrec <- function(prob) {
    poolNrec <-function(mat) {
      n = ncol(mat)-1
      E_f = colSums(t(mat) * 0:n/(2*n))
      n * (1 - 2 * E_f)
    }
    
    nrecs <- vapply(prob, poolNrec, numeric(nrow(prob[[1]])))
    nrecs <- rowSums(nrecs)
  }

  nrecs <- lapply(list(nrec = probasInferred, nrecS = probasSimulated), computeNrec)

  res <- data.table(
    contig = uContigs, d = d[uContigs],
    n = nZ[uContigs, get(poolNames[1]) + ni[2]-get(poolNames[2]) + get(poolNames[3]) + ni[4]-get(poolNames[4])], # n is the true number of recombinants
    do.call(data.table, prop), do.call(data.table, nrecs)
  )

  cat(".")
  res
}


# Evaluating the performance of the inference of nrec and of haplotype phasing ----------------------------

# for this, we simulate contigs at various distances to the SDR.
ds <- 0:20 / 40
res <- rbindlist(lapply(ds, simulate))

# we will only consider contigs with SNPs in both families
twoFams <- res[, !is.na(fam1.prop) & !is.na(fam2.prop)]

# we also retrieve the number of informative SNPs per contig for both families
SNPsPerContig <- lapply(orderedContigs, tabulate)
res[, c("nfam1", "nfam2") := .(SNPsPerContig[[1]][contig], SNPsPerContig[[2]][contig])]

# we make figure S1
library(ggplot2)
library(gridExtra)

temp <- ggplot(data = res[twoFams & n <= 20L, .(n, fam1.prop, fam2.prop, nrec, nrecS, nfam1, nfam2)], mapping = aes(x = factor(n)))

p1 <- temp + geom_abline(intercept = -1, slope = 1, col = "darkgrey") + 
  geom_violin(aes(y = nrec), fill = rgb(1, 0.6, 0.6), col = NA, scale = "width") + xlab("") +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.45), "cm")) + 
  ylab("nrec (phased haplotypes)")

p2 <- temp + geom_abline(intercept = -1, slope = 1, col = "darkgrey") + 
  geom_violin(aes(y = nrecS), fill = rgb(1, 0.6, 0.6), col = NA, scale = "width") + xlab("") +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.45), "cm")) + 
  ylab("nrec (true haplotypes)")

p3 <- temp + geom_violin(
  aes(y = (fam1.prop * nfam1 + fam2.prop * nfam2) / (nfam1 + nfam2), weight = nfam1 + nfam2), 
  fill = rgb(1, 0.6, 0.6), scale = "width", col = NA
  ) +
  xlab("simulated number of recombinants") + ylab("proportion of poperly phased SNPs") +
  theme(plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))

pdf("../../figureS1.pdf", 6, 12)
grid.arrange(p1, p2, p3, nrow = 3)
dev.off()
rm(temp, p1, p2, p3)



# Evaluating the power of the assignments to sex chromosomes ------------------------

# to do so, we simulate autosomal contigs (1000 simulations per contig)
# this is done on a compute server due to the CPU time and RAM required
res <- mclapply(rep(0.5, 1000), simulate, mc.cores = 5)
saveRDS(res, "simuAutosomes_1000repl.RDS")

# we also save nrec only, as we transfer results to a PC with less RAM and over
# a slow internet connection (thanks Covid) this is also why we use the RDS
# format and convert nrec to integer (32 bits instead of 64 for real numbers)
# these results are also used to assign "real" contigs to sex chromosomes (at
# step 05)
res = rbindlist(res)
saveRDS(as.integer(res$nrec * 10^7), stri_c("nrec_1000repl.RDS"))

saveRDS(as.integer(res$nrecS * 10^7), "nrecS_1000repl.RDS")

# we now are on the home PC
# we import all the values of nrec obtained by the simulations assuming autosomes
nrecs <- lapply(list.files(pattern = "nrec[S]*_1000repl"), readRDS)

# we use these results to generate different nrec quantiles for each contig
# these quantiles correspond to the significance levels we investigate

probs <- c(1/1000, 5/1000, 1/100, 5/100)

# we prepare column names corresponding to the above quantiles
qNames <- c("q001", "q005", "q01", "q05")

nrecQuantiles <- function(nrec) {
  # we add a column indicating the contig for each value of nrec
  simu <- data.table(contig = rep(uContigs, 1000), nrec = nrec / 10^7)

  # and we compute the quantiles for each contig
  stats <- simu[, .(q = quantile(nrec, probs)), by = contig]
  stats[, qName := rep(qNames, length.out = .N)]
  
  # we place the different quantiles in different columns
  stats <- dcast(stats, contig ~ qName, value.var = "q")
  stats[match(res$contig, contig), -"contig"]
}

quantiles <- lapply(nrecs, nrecQuantiles)


# We compute the proportion of contigs assigned to sex chromosomes
# depending on d, the significance level an whether SNP type was inferred

# we rbind the data we need for both simulations in a single table
bothRes <- rbind(
  data.table(res[twoFams, .(d, nrec)], quantiles[[1]][twoFams], inferred = T),
  data.table(res[twoFams, .(d, nrec = nrecS)], quantiles[[2]][twoFams], inferred = F)
)

# which facilitates the computation of proportions
props <- bothRes[, .(
  q001 = mean(nrec < q001), q005 = mean(nrec < q005),
  q01 = mean(nrec < q01), q05 = mean(nrec < q05)
), by = .(d, inferred)]

setorder(props, inferred, d)

# we also estimate the proportion of assigned contigs among those at or below a
# certain distance to the SDR
cumu <- lapply(ds, function(x) {
  bothRes[d <= x, .(
    q001 = mean(nrec < q001), q005 = mean(nrec < q005),
    q01 = mean(nrec < q01), q05 = mean(nrec < q05)
  ), by = inferred]
})

cumu <- data.table(d = rep(ds, each = 2), rbindlist(cumu))
setorder(cumu, inferred, d)

pdf("../../figureS3.pdf", 7, 9)

par(mfrow = 2:1, lwd = 1, mai = c(0.8, 0.8, 0.4, 0.4))
ggBackground(range(ds * 100), 0:1, xlab = "distance to the SDR (cM)", ylab = "proportion assigned to sex chromosomes")
with(props[inferred == T, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 1, lwd = 1.5))
with(props[inferred == F, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 3, lwd = 1.5))

legend(
  "topright", lty = c(1, 1, 1, 1, 1, 3), col = c(1:4, 1, 1), 
  legend = c(probs, "inferred SNP type", "true SNP type"), 
  cex = 0.7, bg = grey(0.92), box.col = "white"
  )

text(-7, 1.1, "A", xpd = NA, cex = 1.2, font = 2)

ggBackground(range(ds * 100), 0:1, xlab = "highest distance to the SDR (cM)", 
             ylab = "proportion assigned to sex chromosomes")
with(cumu[inferred == T, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 1, lwd = 1.5))
with(cumu[inferred == F, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 3, lwd = 1.5))
text(-7, 1.1, "B", xpd = NA, cex = 1.2, font = 2)

dev.off()

diffProps <- dcast(cumu[, .(d, inferred, q001)], 
                   d ~ ifelse(inferred, "YES", "NO"), value.var = "q001")
diffProps[which.max(abs(YES - NO))]




# investigating the reduction of crossover rates near the SDR -------------------------------
# we simulate contigs assuming uniform recombination rates
# For this, one simulation is enough (each contig is assigned only one d)
res <- simulate(d = runif(nCont, min = 0, max = 0.5))

# we only save nrec, which we use to make figure 4 at stage 05
saveRDS(res$nrec, "../simu/nrecContinuousCrossovers.RDS")


# same as above, but in case the SDR is at 7 cM from the edge of a chromosome (as suggested by results)
cM = runif(nCont*2, min = 0, max = 0.5)

# we double the probability to sample a contig that is less than 7 cm away from the SDR
probs = rep(1, nCont*2)
probs[cM < 0.07] = 2
probs = probs / sum(probs)
cM = sample(cM, size = nCont, prob = probs)

res <- simulate(d = cM)
saveRDS(res$nrec, "nrecContinuousCrossovers_7cMFromEdge.RDS")


