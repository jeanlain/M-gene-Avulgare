##%######################################################%##
#                                                          #
####           We phase maternal haplotypes             ####
####      and estimate haplotype frequencies using      ####
####               data from the F1 pools               ####
#                                                          #
##%######################################################%##


# NOTE : we denote the "SDR allele to which the maternal allele is linked",
# which is quite a mouthful, simply as "SNP type", and we use integer numbers to
# denote that (for performance reasons). Hence a SNP is of "type 1" if its
# paternal allele is linked to the M allele and is of "type 2" if this allele is linked
# to the m allele.
# IMPORTANT : to avoid using different functions for daughters and sons and to simply
# the code, we revert SNP type in sons (a type 1 SNP in daughters becomes a
# type-2 SNP in their brothers and reciprocally). This allows using the same
# equations in both sexes, and is equivalent to the equations in the Cordaux et al. 2021 
# (this was actually the solution described in the first version of the manuscript
# before the revision, but this made the biological explanations much more
# complex, so we dropped it in the revision). We refer to the reverted SNP type
# as "SNP group" (in both sexes)

# Because of this, in daughters f is the frequency of the m-linked haplotype, not of the
# M-linked! Like the M-linked haplotype in sons, the m allele should be inherited by all
# daughters for a contig that is totally linked to the SDR (f = 0.5). So we are
# particularly interested in p(f = 0.5) for both sexes, which also simplifies the
# code (we have the same expectations in both sexes)

source("WZ_functions.R")

# our input is the table of informative SNPs generated at step 02
  vcf <- fread("../tables/fam2/vcfAndF1countsLikelihoodQ20_A1_angsd.txt")

# and the scans of F1 bases generated at step 03 (imported later)

nSons = 30L
nDaughters = 30L

# we first filter SNPs --------------------------------------------

# assigning types to SNPs = phasing haplotype (see paper)
# 1 means that the maternal allele is linked to the W
# 2 means that it is linked to the Z.
vcf[, type := 2L - ((Cs/ Rs) >= (Cd / Rd))]

# when frequencies are the same between the sexes, we attribute types at random
# (actually, based on whether the SNP position is even, to permit perfect
# replication in case this code is ran several times)
vcf[(Cd / Rd) == (Cs / Rs), type := POS %% 2L + 1L]

vcf[, p := pFgivenCountsSNPs(Cd, Rd, group = 3L-type, ni = nDaughters) * 
      pFgivenCountsSNPs(Cs, Rs, group = type, ni = nSons)]

# We check if allele frequencies are plausible given parental genotypes. In
# particular, we check if allele 1 is not too frequent in F1s (should be AT MOST
# 50%, considering that a frequency of 0.5 is highly unlikely to being with)
# We do it by computing the p-value of an exact binomial test that f > 0.5
vcf[, plausible := pbinom(Cd - 1L, Rd, 0.5, lower.tail = F) * pbinom(Cs - 1L, Rs, 0.5, lower.tail = F)]

# we flag the SNPs that we retain.
vcf[, retained := Cd + Cs > 0L & Rd + Rs > 0L # those covered in both pools and where the maternal allele was found
& (plausible > 0.05 | pmin(Cd / Rd, Cs / Rs) < 1 / (2*nSons)) &
  Rd + Rs <= quantile(Rd + Rs, 0.95) & # those of sequencing depth not higher than the 95% quantile
  mother.DP >= 15L &
  Rd + Rs > (dA + dC + dG + dT + sA + sC + sG + sT) * 0.75] # the maternal + paternal allele must constitute at lest 75% of the F1 reads

# but we retain those that yield not too bad probability of f = 0.5
vcf[p > 0.05, retained := T]

vcf[is.na(retained), retained := F]


# We locate "outlier" SNPs that give a much lower probability that f=0.5
# than the other SNPs of the same contig.
# remember that f should equal 0.5 even in daugthers if the contig did not recombine 
# with the SDR (see introductory NOTE)
# we reframe the table so as all pools are on different rows and to obtain a
# single column for C and R, but we keep information on the sex with a new
# logical "sex" column which is TRUE for females. Doing so simplifies the script
counts <- vcf[retained == T, data.table(
  CHROM = rep(CHROM, 2), genomePos = rep(genomePos, 2),
  A1 = rep(A1, 2), A2 = rep(A2, 2),
  C = c(Cd, Cs), R = c(Rd, Rs), group = c(3L - type, type),
  female = rep(c(T, F), each = .N)
)]


# so we compute the (-log) posterior probability that f = 0.5 for each SNP
counts[female == T, p := -log(pFgivenCountsSNPs(C, R, group, ni = nDaughters))]
counts[female == F, p := -log(pFgivenCountsSNPs(C, R, group, ni = nSons))]


# and we apply our criterion to locate outlier SNPs
limit <- counts[, pLimit(p), by = .(cr = CHROM, fem = female)]
counts[, lim := limit[chmatch(stri_c(CHROM, female), stri_c(cr, fem)), V1]]

# we discard outliers in the vcf table. If a SNPs is an outlier in a pool, we
# discard it from both pools of a family, but not in the other family (if it is
# present)

vcf[genomePos %in% counts[p > lim, genomePos], retained := F]


setorder(vcf, genomePos)

# we now count c1, r1, c2, r2 and compute P(f|cri) -----------------------------------
# Here, 1 refers to W in daughters and 2 to Z in sons, 
# 2 refers to Z in daughters and to W in sons 
# We do not use the notation cw, cz, rw and rz, as our implementation is not 
# a direct transcription of the method section (see the introductory NOTE)

# we import contig lengths to assign absolute positions to contigs in these tables
contigLengths <- fread("../contigLengths.txt", header = T)

# we make tables of the retained SNPs that we will use for reach pool
# at first, one per family
SNPs = vcf[retained == T, .(CHROM = chmatch(CHROM, contigLengths$contig), 
                                  genomePos, A1, A2, group = type)]

# for daughters, we revert SNP types (see the NOTE at the begining)
SNPd = copy(SNPs)
SNPd[,group := 3L-group]

# we duplicate each table, to have one per pool (2 pools per family)
SNPs = list(SNPd, SNPs)



# we import the scans of F1 bams, which is a list of data.tables
scans <- list.files("../BFog3054/1-F1", pattern = "basesAtAngsdSNPs.rds", full.names = T)
F1scans <- lapply(scans, readRDS)

# we name the tables
names(F1scans) <- ifelse(grepl("female|daug", scans), "daughters", "sons")
F1scans = F1scans[c("daughters","sons")]

# we obtain the cri variables (read counts for different SNP groups)
# see the countReads() function in WZ_functions.R
readCounts <- setNames(
  Map(countReads, F1scans, SNPs, minMapq = 20L, minPhred = 20L, nrow(contigLengths), T, names(F1scans)), 
  names(F1scans)
  )

# we reclaim some RAM
rm(F1scans, counts, SNPs) ; gc()


# we compute the posterior probabilities of all possible haplotype frequencies
probs <- Map(function(counts, ni) do.call(pFgivenCountsHapl, c(counts, ni)), readCounts, c(nDaughters, nSons))


# we compute nrec :
# First, the expected value of f (weighted mean of f values given their posterior)

E_f <- vapply(probs, function(mat) colSums(t(mat) * (1:ncol(mat)-1)/(2*(ncol(mat)-1))), numeric(nrow(probs[[1]])))



# then nrec
nrec <- nDaughters * (1 - 2 * E_f[,"daughters"]) + nSons * (1 - 2 * E_f[,"sons"])

# we also record the highest posterior probability for f in each contig
maxP <- lapply(probs, rowMaxs)

# and the probability that f = 0.5 for each contig
p05 <- lapply(probs, function(x) x[, ncol(x)])

# we add these to the read count tables
readCounts <- Map(data.table, readCounts, pM = maxP, p05 = p05)

# and we combine the results for all pools
probs <- data.table(
  contigLengths,
  do.call(data.table, readCounts),
  nrec,
  p = Reduce("*", p05)
)

writeT(probs, "../tables/fam2/countsAndPosteriors.txt")
writeT(vcf, "tables/vcfAndF1countsLikelihoodQ20_A1_angsd.txt")

writeT(probs[p > 0.5, .(contig, 1, length)], "../tables/Mcandidates.bed", col.names = F)

probs2 = fread("/Users/jean/Documents/scrap/Mm/BFog3054/tables/countsAndPosteriors.txt")

probs = merge(probs, probs2, by = c("contig", "length"), all = T, suffixes = c(".fam1", ".fam2"))

probs[,c("nrec","p") := .(nrec.fam1 + nrec.fam2, p.fam1 * p.fam2)]


# checking that individuals from other lineages lack the paternal allele of type 1-SNPs.----------------------------------
# if they don't, it means that these SNPs cannot be responsible for the male phenotype

setwd("~/jpeccoud_data/genomeAV")
bams = list.files(pattern = ".bam$", recursive = T)
bams = bams[!grepl("BH|Mm|F1", bams) & file.size(bams) > 10e9]

res = mclapply(bams, function(bam) {
  system(paste("samtools view -h -b -ML ~/jpeccoud_data/genomeAV/Mm_SNPs/Mcandidates.bed", bam, ">", gsub(".bam$",".Mcandidates.bam", bam)))
}, mc.cores = 4)


sel = vcf[, retained ==T & type == 1L & CHROM %chin% probs[p > 0.5, contig]]

# we save a bed file of the type-1 SNPs
bed = vcf[sel, .(CHROM, POS)]
writeT(bed, "/Users/jean/Documents/scrap/Mm/others/type1SNPs.bed", col.names = F)

setwd("/Users/jean/Documents/scrap/Mm/others")

# bam files from other individuals, only covering the selected contigs (otherwise the pileup would be very slow)
bams = list.files(pattern = ".bam$")

# we generate pileups to count bases
# somehow, a single pileup would only return the first 243 positions of the bed
test = mclapply(bams, function(bam) system(stri_c("samtools mpileup -q 20 -a -B -l type1SNPs.bed ", bam, " > ", bam, ".pileup.txt")), mc.cores = 4)

pile = lapply(list.files(pattern = "pileup.txt"), fread, header = F, select = 5)
pile = do.call(data.table, c(bed, pile))
setnames(pile, 3:ncol(pile), c("mWXA","mZM","fWXA","fZM","mP34","dBF","sBF","mWXf"))

# we count the bases at each position
regs = setNames(c("a|A", "c|C","g|G", "t|T"),c("A","C","T","G"))

counts = lapply(pile[,-(1:2), with = F], function(cigar) {
    vapply(regs, function(reg) stri_count(cigar, regex = reg), integer(length(cigar)))
  }
)

# we compute the sequencing depth
dep = sapply(counts, rowSums)
A1 = vcf[retained ==T & type == 1L & CHROM %chin% probs[p > 0.5, contig], A1]
nr = length(A1)
nA1 = vapply(counts, function(count) count[cbind(1:nr, A1)], integer(nr))

res = data.table(vcf[sel, .(CHROM, POS, father.AD, mother.AD, A1, A2, Cd, Rd, Cs, Rs)], nA1)

ok = res[rowMaxs(nA1 / dep) < 0.05]



contigs = ok[,table(CHROM)]
contigs[sort(names(contigs))]

test = system(paste("samtools mpileup -q 20 -a -B -l type1SNPs.bed", paste(bams, collapse = " "), "> type1SNPs.pileup.txt"))



