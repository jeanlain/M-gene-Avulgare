
##%######################################################%##
#                                                          #
####    We locate genomic regions which are suitable    ####
####   for the SDR, taking advantage of both families   ####
#                                                          #
##%######################################################%##

# we reconstruct parental haplotypes to find SNPS that have recombined
# with the SDR during the divergence of the WXA and ZM lines
# to do this, we retrieve useful SNPs in the parents

setwd("/Users/jean/Library/CloudStorage/OneDrive-UniversitédePoitiers/recherche/sex_determinism/symchrosex/M_gene/scripts_SNPs")
source("WZ_functions.R")

# we import results from earlier steps, as we will select SNPs on sex chromosomes (otherwise, there would be too many SNPs)
probs = fread("../tables/countsAndPosteriors.txt")

# we import position of SNPs, to select those in the contigs that are on the sex chromosomes

pos = fread("../0-F0/fam1_2_mother_father_remapped_gatk.pos", header = T, drop = 3)

# because of the large number of SNPs, we only retained those on contigs assigned to the M chromosome
retained = pos$chr %chin% probs[WZ == T, contig]
pos = pos[retained == T,]


# we import the beagle
beagle = fread("../0-F0/fam1_2_mother_father_remapped_gatk.beagle", header = T, drop = 1)
beagle = beagle[retained == T,]
newNames =   paste(rep(c("Mo_", "Fa_"), each = 3), rep(c("fam1", "fam2"), each = 6),  rep(c("_Maj","_Het", "_Min"), 4), sep = "")

setnames(beagle, 3:ncol(beagle), newNames)


# we import the read counts
counts = fread("../0-F0/fam1_2_mother_father_remapped_gatk.counts", header = T)
counts = counts[retained == T,]
counts[,V17 := NULL]
newNames =   paste(rep(c("Mo_", "Fa_"), each = 4), rep(c("fam1", "fam2"), each = 8),  c("_A","_C", "_G","_T"), sep = "")
setnames(counts, newNames)

vcf = data.table(pos, beagle, counts)
setnames(vcf, 1:4, c("CHROM","POS","A1", "A2"))
vcf[,A1 := A1+1L]
vcf[,A2 := A2+1L]




genotypes = function(ind, fam) {
  col = grep( paste(ind, fam, sep = "_"), names(vcf), fixed = T)
  vcfOneInd = vcf[,c(3:4, col), with = F]
  setnames(vcfOneInd, 3:ncol(vcfOneInd), c("Maj","Het","Min", "A","C","G","T"))
   
  scores = vcfOneInd[,cbind(Maj, Min, Het)]
  maxScores = rowMaxs(scores)
  geno = rowMatches(maxScores, scores)  # genotype: 1 will mean A1/A1, 2 = A2/A2,  3 = A1/A2
  vcfOneInd[geno == 1L, c("allele1","allele2") := .(A1, A1)]
  vcfOneInd[geno == 2L, c("allele1","allele2") := .(A2, A2)]
  vcfOneInd[geno == 3L, c("allele1","allele2") := .(A1, A2)]
  bases = vcfOneInd[,cbind(A,C,G,T)]
  rows = 1:nrow(bases)
  countsA1 = bases[cbind(rows, vcfOneInd$A1)]
  countsA2 = bases[cbind(rows, vcfOneInd$A2)]
  
  vcfOneInd[, cbind(maxScores, geno, allele1, allele2, countsA1, countsA2)]
  
}

alleles = Map(genotypes, ind = rep(c("Mo","Fa"), 2), fam = rep(c("fam1","fam2"), each = 2))
alleles = do.call(cbind, alleles)
maxScoreCols = colnames(alleles)=="maxScores"
allelesCols = colnames(alleles) %chin% c("allele1","allele2")

nHaplA1 = rowSums(alleles[,allelesCols] == vcf$A1, na.rm = T)
nHaplA2 = rowSums(alleles[,allelesCols] == vcf$A2, na.rm = T)
minMaxScores = rowMins(alleles[,maxScoreCols])

vcfRetained = data.table(vcf[,.(CHROM, POS, A1, A2)], alleles[,!maxScoreCols & !allelesCols])

newNames = paste(rep(c("Mo", "Fa"), each = 3), rep(c("fam1", "fam2"), each = 6), colnames(vcfRetained)[5:ncol(vcfRetained)], sep ="_")
setnames(vcfRetained, 5:ncol(vcfRetained), newNames)
vcfRetained[,keep := nHaplA1 >= 2L & nHaplA2 >= 2L & minMaxScores >= 0.9]


contigLengths = fread("../contigLengths.txt")
contigStarts <- setNames(c(0L, cumsum(contigLengths$length[-nrow(contigLengths)])), contigLengths$contig)
vcfRetained[, genomePos := POS + contigStarts[CHROM]]


informativeSNPSFiles = list.files("../tables", pattern = "vcfAndF1countsLikelihoodQ20_A1_angsd.txt", recursive = T, full.names = T)

informativeSNPs = lapply(informativeSNPSFiles, fread, header = T, 
                         select = c("CHROM", "POS", "genomePos", "A2","A1","A2_mother_count","A1_mother_count","A2_father_count","A1_father_count", "type", "retained"), 
                         col.names = c("CHROM", "POS", "genomePos", "A1","A2","Mo_fam1_countsA1","Mo_fam1_countsA2","Fa_fam1_countsA1","Fa_fam1_countsA2", "type", "retained"))


reshapeVCF = function(vcf) {
  vcf[, c("Mo_fam1_geno","Fa_fam1_geno") := .(1L, 3L)]
  vcf = vcf[CHROM %chin% probs[WZ == TRUE, contig]]
}

informativeSNPs = lapply(informativeSNPs, reshapeVCF)

setnames(informativeSNPs[[2]], gsub("fam1","fam2", colnames(informativeSNPs[[2]])))

merged = merge(informativeSNPs[[1]], informativeSNPs[[2]], by = c("CHROM", "POS", "genomePos"), all = TRUE, suffixes = c("_fam1","_fam2"))

vcfRetained = merge(vcfRetained, merged[,.(genomePos, type_fam1, type_fam2, retained_fam1, retained_fam2)], by = "genomePos", all.x = T)

merged = merged[!genomePos %in% vcfRetained$genomePos]

alleles = merged[,cbind(A1_fam1, A2_fam1, A1_fam2, A2_fam2)]
nTimes = sapply(1:4, function(base) rowSums(alleles == base, na.rm = T))
nAlleles = rowSums(nTimes  > 0L)
merged = merged[nAlleles < 3L,]
pb = merged[, A1_fam1 != A1_fam2]
merged[pb, c("A1_fam2","A2_fam2","Mo_fam2_countsA1","Mo_fam2_countsA2","Fa_fam2_countsA1","Fa_fam2_countsA2", "Mo_fam2_geno") := 
         .(A2_fam2, A1_fam2, Mo_fam2_countsA2, Mo_fam2_countsA1, Fa_fam2_countsA2, Fa_fam2_countsA1, 2L) ]

merged[,A1 := ifelse(is.na(A1_fam1), A1_fam2, A1_fam1)]
merged[,A2 := ifelse(is.na(A2_fam1), A2_fam2, A2_fam1)]

merged[,c("A1_fam1","A2_fam1", "A1_fam2","A2_fam2") := NULL ]

vcfRBind = rbind(vcfRetained, data.table(merged, keep = T))
vcfFinal = vcfRBind[keep == T | !is.na(retained_fam1) | !is.na(retained_fam2), -"keep"]

posFile = "/Users/jean/Documents/scrap/Mm/bothFam/retainedSNPs.txt"
writeT(vcfFinal, posFile)


# we retrieve the bases in descendants, as we will check the reliability of SNPs -------------------
library(parallel)
setwd("~/Documents/scrap/Mm")
bamFiles = list.files("1-F1", "genome_markdup.bam$", full.names = T, recursive = T)
posFile = list.files("0-F0", "retainedSNPs.txt$", full.names = T, recursive = T)

# we launch the script to scan each F1 bam in parallel (using 5 CPUs for each)
m <- mcMap(
  function(posFile, bamFile) {
    system(paste("Rscript scanBamAtPositionsAngsd.R", posFile, bamFile, 20))
  },
  rep(posFile, each = length(bamFiles)), # each position file is used for two bam files (two pools)
  bamFiles,
  mc.cores = length(bamFiles)
)


# we list the resulting files
scans <- list.files("~/Documents/scrap/Mm", pattern = "basesSexChromsomesSNPs.rds$", full.names = T, recursive = T)

# we import the scans of F1 bams, which are data.tables
F1scans <- lapply(scans, readRDS)
names(F1scans) <- paste(ifelse(grepl("females|daughters", scans), "daughters", "sons"), rep(c("fam1","fam2"), each = 2), sep = "_")

setorder(vcfFinal, genomePos)
SNPpos <- vcfFinal$genomePos

# and we count reads
f1Counts <- Map(baseMatrix, scannedBam = F1scans, genomePos = list(SNPpos, SNPpos, SNPpos, SNPpos), minMapQ = 20L, onlyProper = F)

veriSNPs = function(fami) {  # verifies the reliability of SNPs based on read counts from parents and F1s
  baseMatrices = f1Counts[grep(fami, names(f1Counts))]
  
  # we first flag all SNPs at position with excessive sequencing depth in each parent or pool
  depths = sapply(baseMatrices, rowSums)
  cols = which(grepl(fami, colnames(vcfFinal)) & grepl("counts", colnames(vcfFinal)))
  parentCounts = vcfFinal[, as.matrix(cols), with = F]
  depths = data.table(depths, parentCounts[,1] + parentCounts[,2], parentCounts[,3] + parentCounts[,4])
  quantiles = sapply(depths, quantile, probs = 0.95, na.rm = T)
  excessDepth = t(t(as.matrix(depths)) > quantiles)
  excessDepth = rowSums(excessDepth, na.rm = T) >= 1L

  F1counts = Reduce("+", baseMatrices)  # from now on, we pool F1s with respect to read counts
  depth = rowSums(F1counts)           # the sequencing depth other all F1s
  
  cols = which(grepl(fami, colnames(vcfFinal)) & grepl("geno", colnames(vcfFinal)))
  vcfFam = vcfFinal[, cols, with = F]
  setnames(vcfFam, c("Mo_geno","Fa_geno"))
  vcfFam = data.table(vcfFinal[,.(A1, A2)], vcfFam)
  
  # checks if a base for which at least one parent is homozygous is frequent enough in F1s
  # this base should have a frequency of at least 0.5 in F1s.
  A1fixed = vcfFam[, Mo_geno == 1L | Fa_geno == 1L]
  A2fixed = vcfFam[, Mo_geno == 2L | Fa_geno == 2L]
  A1fixed[is.na(A1fixed)] = F
  A2fixed[is.na(A2fixed)] = F
  
  dominantBase = vcfFam$A1  
  dominantBase[A2fixed] = vcfFam$A2[A2fixed];
  rows = 1:nrow(F1counts)
  dominantCounts = F1counts[cbind(rows, dominantBase)]
  
  p = pbinom(dominantCounts, depth, 0.5)
  tooRare = (A1fixed | A2fixed) & p < 0.05
  
  # checks that an allele that is fixed in parents is almost exclusively found in F1
  
  baseFreq = dominantCounts / depth
  fixed = vcfFam[,(Mo_geno == 1L & Fa_geno == 1L) | (Mo_geno == 2L & Fa_geno == 2L) ]
  fixed[is.na(fixed)] = F
  fixedBaseToRare = fixed & (baseFreq < 0.95 | depth < 15L | depth - dominantCounts > 3L)
  
  # checks if an allele for which parent are homozygous (and different) has a frequency close to 0.5 in F1s
  
  bilateralBinomTest = function(x, n, p = 0.5) {
    x = pmin(x, n-x)
    pbinom(x, n, 0.5)*2
  }
  
  p = bilateralBinomTest(dominantCounts, depth)
  
  fixedParent = A1fixed & A2fixed
  fixedParentProblem = fixedParent & (p < 0.01)
  
   
  # check heterozygous genotypes in parents are consistent with read counts
  # (either base should not be too rare in parent reads)
  
  p = Map(bilateralBinomTest, parentCounts[,c(1, 3), with = F], depths[,3:4, with = F])
  p = do.call(cbind, p)
  
  hetero = vcfFam[, cbind(Mo_geno == 3L, Fa_geno == 3L)]
  heteroParentProblem =  rowSums(p < 0.01 & hetero) > 0L
  
  cbind(excessDepth, tooRare, fixedBaseToRare, fixedParentProblem, heteroParentProblem)
  
}

checks = lapply(c("fam1","fam2"), veriSNPs)
famCheck = lapply(checks, function(x) rowSums(x, na.rm = T) == 0L)
vcfFinal[, c("ok_fam1","ok_fam2") := famCheck]


# To select the relevant contigs, we import our main table of results ---------------
setwd("/Users/jean/Library/CloudStorage/OneDrive-UniversitédePoitiers/recherche/sex_determinism/symchrosex/M_gene/scripts_SNPs")
probs = fread("../tables/countsAndPosteriors.txt")

sdr = vcfFinal[CHROM %in% probs[p > 0.5, contig]]

# we select SNPs according to our checks. We could have done that earlier
# but is is better to filter late in case we need to go back
# we exclude SNPs that were considered unreliable at step 03
sdr = sdr[!(!is.na(retained_fam1) & retained_fam1 == F) &
            !(!is.na(retained_fam2) & retained_fam2 == F)]

# we exclude SNPs that did not pass our checks,
# but we don't exclude informative SNPs that we retained at step 03
# for non-informative SNPs, fathers must be homozygous
sdr = sdr[(retained_fam1 == T | (ok_fam1 & Fa_fam1_geno != 3L)) & 
              (retained_fam2 == T | (ok_fam2 & Fa_fam2_geno != 3L)) ]  

# we also exclude rare SNPs that did not produce consistent results between allele calls
# i.e., informative SNPs for which the father was not found heterozygous this time
# or the mother not found homozygous
pb = sdr[, retained_fam1 == T & (Fa_fam1_geno != 3L | Mo_fam1_geno == 3L)]
pb[is.na(pb)]= F
pb2 = sdr[, retained_fam2 == T & (Fa_fam2_geno != 3L | Mo_fam2_geno == 3L)]
pb2[is.na(pb2)] = F
sdr = sdr[!pb & !pb2,]


# to infer recombinant and causal SNPs, we use data from other lineage that don't possess the M allele
# as we will use mpileup to count bases at these SNPs in various bam files, 
# we first extract data from the relevant contig, which we save in a bed
writeT(probs[p > 0.5, .(contig, 1, length)], "../tables/Mcandidates.bed", col.names = F)

setwd("~/jpeccoud_data/genomeAV")
bams = list.files(pattern = ".bam$", recursive = T)
bams = bams[!grepl("BH|Mm|F1|assembly", bams) & file.size(bams) > 10e9]
outdir = "~/jpeccoud_data/genomeAV/Mm_SNPs/Mcandidate_other_lineages"
bed = paste(outdir, "Mcandidates.bed", sep ="/")
if(!dir.exists(outdir)) {
  dir.create(outdir)
}

out =  paste(outdir, gsub(".bam$",".Mcandidates.bam", basename(bams)), sep = "/")

res = mcMap(function(bam, out) {
  system(paste("samtools view -h -b -ML", bed, bam, ">", out))
},bams, out, mc.cores = 4)


# we save a bed file of SNPs from these contigs
bed = sdr[, .(CHROM, POS)]
outdir = "~/Documents/scrap/Mm/Mcandidate_other_lineages"
bedFile = paste(outdir, "candidateContigsSNPs.bed", sep ="/")
writeT(bed, bedFile, col.names = F)

bams = list.files(outdir, pattern = ".bam$", full.names = T)

# we generate pileups to count bases
pile = mclapply(bams, function(bam) system(paste("samtools mpileup -q 20 -a -B -l ", bedFile, " ", bam, " > ", bam, ".pileup.txt",  sep = "")), mc.cores = 4)


pileups = list.files(outdir, pattern = "pileup.txt", full.names = T)
pile = lapply(pileups, fread, header = F, select = 5)

pile = do.call(data.table, c(bed, pile))
setnames(pile, 3:ncol(pile), gsub(".Mcandidates.bam.pileup.txt","", basename(pileups)))


# we count the bases at each position
regs = setNames(c("a|A", "c|C","g|G", "t|T"),c("A","C","G","T"))

counts = lapply(pile[,-(1:2), with = F], function(cigar) {
  vapply(regs, function(reg) stri_count(cigar, regex = reg), integer(length(cigar)))
})


# we compute the sequencing depth
dep = sapply(counts, rowSums)

# we counts reads for alleles A1 and A2
rows = 1:nrow(sdr)
nA1 = vapply(counts, function(count) count[cbind(rows, sdr$A1)], integer(length(rows)))
nA2 = vapply(counts, function(count) count[cbind(rows, sdr$A2)], integer(length(rows)))

A1Present = ((pbinom(nA1, dep, 0.5) > 0.05)  & nA1 >= 2L) * 1L
A2Present = ((pbinom(nA2, dep, 0.5) > 0.05) & nA2 >= 2L) * 2L

geno = A1Present + A2Present # geno of value 1 == homozygous for A1, 
# 2 == homozygous for A2, and 3 == heterozygous (0 == undetermined)

geno[dep < 5 | (nA1 + nA2 < dep*0.9)] = 0L
geno = data.table(sdr[, .(Fa_fam1_geno, Mo_fam1_geno, Fa_fam2_geno, Mo_fam2_geno)], as.data.table(geno))

carriedAlleles = function(genotype) {
  alleles = sdr[, cbind(a1 = A1, a2 = A1)]
  alleles[genotype == 2L, ] = sdr[genotype == 2L, cbind(A2, A2)]
  alleles[genotype == 3L, ] = sdr[genotype == 3L, cbind(A1, A2)]
  alleles[genotype == 0L,] = NA
  as.data.table(alleles)
}

alleles = lapply(geno, carriedAlleles)

# for informative SNPs, we make sure that the paternal allele on the left column ("a1")
# is the one linked to the M allele (the allele of the 2nd column would be the one linked to m)
# we use SNP type to infer that
phaseFather = function(dt, type, motherGenotype) {
  dt = copy(dt)
  # if SNP is of type 1, the M-linked allele is the one not carried by the mother
  # if the mother is of genotype 1 (it has allele A1), we must thus take allele A2
  # if the mother is of genotype 2, we take allele A1
  # in fact, if the genotype of the mother corresponds to SNP type, a1 (M-linked) is A2, and a2 is A1
  # for other SNPs, a1 is already A1 and a2 is already A2
  
  dt[type == motherGenotype, c("a1", "a2") := .(a2, a1)]
  dt 
}

alleles$Fa_fam1_geno = phaseFather(alleles$Fa_fam1_geno, sdr$type_fam1, sdr$Mo_fam1_geno)
alleles$Fa_fam2_geno = phaseFather(alleles$Fa_fam2_geno, sdr$type_fam2, sdr$Mo_fam2_geno)

alleleMat = as.matrix(do.call(data.table, alleles))

# we determine if SNPs have recombined with the SDR. To do so, we create 2-SNP
# haplotypes where allele 1 = M and 2 = m
# so the haplotype is encoded as a 2-digit integer (one digit per allele)
MmHapl = t(c(1, 2, 2, 2, 1, 2, 2, 2, rep(2, ncol(dep)*2))*10 + 
             t(alleleMat))

nHapl = countHapl(MmHapl)
# we add an integer column indicating the category of SNPs
# 1 (TRUE) indicates a recombinant SNP
sdr[,rec := as.integer(nHapl == 4L)] 

# 2 when a SNP shows no evidence for recombination and is not informative in both families
sdr[rec != 1L & (is.na(retained_fam1) | is.na(retained_fam2)) , rec := 2L]       

# 3 if the SNP is informative in both families
sdr[rec != 1L & !is.na(retained_fam1) & !is.na(retained_fam2), rec := 3L]       

# 4 for potentially causal SNPS
sdr[nHapl == 2L & type_fam1 == 1L & type_fam2 == 1L, rec := 4L]    


# we save the SNP category to disk as we will need it for figure 7
writeT(sdr[, .(CHROM, POS, rec)], "tables/SNPcategory.txt")


# we assign SNPs to recombination blocks -----------------------------------------

dt = data.table(sdr[,.(rec, CHROM, genomePos)], do.call(data.table, alleles[1:4]))
setorder(dt, CHROM, genomePos)

# we retrieve coordinates of genomic blocks (which are those of SNPs at the edges)
breaks = dt[, breakPoints( data.table(CHROM, .SD)), by = CHROM]$V1


contigLengths = readNamedVector("../contigLengths.txt")

# we sort the coordinates in ascending order to assign SNPs to blocks
# we do this by creating bins (not forgetting to add the two extreme coordinates)
breaks = sort(unique(c(1L, breaks, sum(contigLengths)+1L)))
sdr[,block := .bincode(genomePos, breaks = breaks, right = F)]

# in case a bin overlaps two contigs, we make sure it does not correspond to a unique block
contigBlock = sdr[,paste(CHROM, block)]
sdr[,block:= chmatch(contigBlock, unique(contigBlock))]


# we get relevant metrics for blocks
blocks = sdr[,data.table(
    start = min(POS), end = max(POS),              # positions of first and last SNP 
    data.table(rbind(tabulate(rec, nbins = 4))),   # and number of SNPs per category
    used = sum(!is.na(rec))                         # and number of SNPs that inform on recombination
    ), by = .(CHROM, block)]

# we select blocks composed of more than one SNP or that contain causal SNPs
blocks = blocks[start < end | V4 > 0L]


# we retrieve contigs that did not have SNPs that could be used to assess
# recombination with the SDR (hence not in the genotypes table). We will consider
# that there is no proof of recombination with the SDR for these contigs so we
# add them as single blocks constituting whole contigs, with just one SNP of rec
# category 2 (non-shared informative SNP not indicating recombination)
missingContigs = setdiff(probs[p > 0.5, contig], blocks$CHROM)

if (length(missingContigs)> 0L) {
  blocks = rbind(blocks, data.table(
    CHROM = missingContigs, block = 1L, start = 1L, end = 1L, 
    V1 = 0L, V2 = 1L, V3 = 0L, V4 = 0L, used = 1L
    ))
}

# we compute block start and end positions (as blocks must be contiguous)
n = 2:nrow(blocks)

# block boundaries are the midpoint between the SNP coordinates of adjacent blocks
blocks[,end := c(as.integer((end[n-1] + start[n])/2), end[.N])]
blocks[, start := c(1L, end[n-1]+1L)]

# the start of the first block of a contig is moved to 1 and the end of the last
# block is moved to the length of a contig
blocks[!duplicated(CHROM), start := 1L]
blocks[!duplicated(CHROM, fromLast = T), end := contigLengths[CHROM]]


# we categorize blocks. We attribute number 1 whenever a block harbors 2 or more
# than 50% recombinant SNPs, otherwise we attribute:
# 2 if it harbors no shared informative SNPs
# 3 if it harbors shared non-causal SNPs
# 4 if it harbors causal SNPs
blocks[, cat := ifelse(V1 >= 2L | V1 / used >= 0.5, 1L, ifelse(V4 > 0L, 4L, ifelse(V3 > 0L, 3L, 2L)))]
blocks[is.na(cat), cat := 2L]

# the proportion of blocks harboring shared informative SNPs
# among recombinant blocks
blocks[,mean(V1 > 1 & V3 + V4 > 0)]

# as we will plot longer contigs first, we retrieve contig lengths
blocks[,len := contigLengths[CHROM]]
setorder(blocks, -len, CHROM)

# blocks will be plotted as rectangles, their top and bottom
# edges are simply their start and end positions in the contigs 
#(-1 for start to ensure that rectangles are contiguous)
blocks[,c("bottom", "top") := .((start-1L)/1000, end/1000)]

# the X position of the left edges are determined by their 
# host contig (+1 per new contig)
row = 2:nrow(blocks)
blocks[,left := cumsum(c(1L, CHROM[row] != CHROM[row-1]))]

# we prepare the second barplot (the inset) showing total lengths
totLengths = blocks[cat != 0L,.(len = sum(end-start+1L)/1000), by = cat]
setorder(totLengths, cat)

# the explanatory labels of this barplot
labels = c("contains SNPs indicating recombination with the SDR",
           "contains no SNPs that are informative in both families",
           "contains SNPs that are informative in both families \nbut no causal SNP",
           "contains potentially causal SNPs")

# this second barplot will appear in the white space of the main plot
# between these Y coordinates
maxY = 250
minY = 150

# so we need to compute the top and bottom of sectors so they fit in there
ratio = (maxY - minY)/sum(totLengths$len)
totLengths[,c("bottom","top") := .(cumsum(c(minY, ratio*len[-.N])), cumsum(ratio*len) + minY)]

# we compute the mid Y position of sectors as we will indicate sizes there
# for the shorter sector, we use the top instead of the middle (as it is too narrow)
totLengths[,mid := ifelse(top -bottom >5, (top + bottom)/2, top)]

# we compute the Y position of the explanatory labels (that are equidistant)
totLengths[,textY := seq(min(mid), max(mid), length.out = .N)]

# the colors corresponding to block categories, from 1 to 4
cols =  c("salmon", "grey", "mediumseagreen","olivedrab2")

pdf("../figure5b.pdf", width = 8, height =7)

par(mai = c(0.2, 0.8, 0.4, 0.35))

# the plot showing all contigs with colors for blocks of different categories
# we prepawre an empty plot to add rectangles to it
with(blocks, plot(
  x = range(left + 1L), y = range(c(top, bottom)), bty = "n", type = "n",
  ylab = "Position in contig (kbp)", las =1,
  xlab = "",  xlim =c(2, max(left) +1),
  xaxt = "n"
))

with(blocks, rect(
  xleft = left, ybottom = bottom, xright = left + 0.9, # we use a 0.1 space between bars
  ytop = top, border = NA, col = cols[cat]
))

# we determine the left X position of the labels of the inset
# so that the longest ends at the X position of the last contig
textR = max(blocks$left) +1L - max(strwidth(labels, cex = 0.7))

# so we can draw the barplot of the inset 
# (we use rectangles as it is more flexible than adding a regular barplot)
xUnit = max(blocks$left)/6
startBar = textR - 2*xUnit
endBar = textR - xUnit
with(totLengths, rect(
  xleft = startBar, ybottom = bottom, xright = endBar, 
  ytop = top, border =NA, col = cols[cat])
  )

# we draw lines connecting the sectors to the captions
l = with(totLengths, Map(function(middle, Y, col) lines(
  x =  endBar + c(0.03, 0.2, 0.8, 1) * xUnit/2, 
  y = c(middle, middle, Y, Y), 
  lend = 1, col = col
  ), mid, textY, cols[cat]) 
)

# we write the explanatory captions
totLengths[, text(x = endBar + xUnit/2, y = textY, labels = labels, pos = 4, col = grey(0.3), cex = 0.7)]

# and the lengths within sectors
totLengths[, text(x = (startBar + endBar)/2, y = mid, labels = paste(round(len), "kbp"), col = grey(0.3), cex = 0.7)]

dev.off()

# checking genes involved
gff = fread("grep -v '^#' ../Annotations.goodScaffIDs.wIntrons.newIDs.gff3", select = c(1, 3:5, 7, 9), header = F,
            col.names = c("Contig","type","Start","End","Strand", "annot"))
gff[, ID := splitToColumns(annot, "_",2)]

genes = gff[type == "gene", -"type"]

library(readxl)
predictions = data.table(read_xlsx("../mGenes_protein_function.xlsx"))
predictions[, id := splitToColumns(`Prot ID`, "_",2)]
genes[, Product := predictions[match(ID, id), `Interpro domains - summary`]]

writeT(genes[CHROM %in% sdr$CHROM], "../tables/genesInSDR.txt")


# we add a column indicating the overlap between each gene and the block with causal SNPs
genes[,ov := intersectsRanges(data.table(CHROM, start, end), blocks[cat == 4L, .(CHROM, start, end)], T)]
writeT(genes[ov > 0L,], "../tables/genesIntersectingBlocks.txt")

# we identify SNPs in CDS
cds = gff[type =="CDS", .(Contig, Start, End, Strand, ID)]
library(Biostrings)
genome = readDNAStringSet("/Users/jean/Documents/scrap/genome.fasta")
bases = c("A","C","G","T") 

coding = sdr[rec == 4L, codingChanges(scaff = CHROM, pos = POS, gff = cds, scaffSeqs = genome, base = bases[A1], base2 = bases[A2])]
coding[, SNP := stri_c(posInGene,": ", refCodon, "->", altCodon, " ", refAA, "->", altAA)]
perGene = coding[refAA != altAA,.(stri_flatten(SNP, collapse = "\n")), by = gene]
genes[,`Causal SNPs` := perGene[match(ID, gene), V1]]
genes[Product == "NA", Product := "Hypothetical protein"]

library(xlsx)
write.xlsx(genes[,-c("annot", "ID")], "../genesInSDR2.xlsx", showNA = F, row.names = F)

writeT(coding, "../tables/causalSNPsInCDS.txt")

writeT(sdr, "../tables/SNPsInCandidateContigs.txt")


# alternative figure showing the location of causal SNPs and of genes on contigs
lengths = contigLengths[contig %in% sdr$CHROM]
setorder(lengths, -length)
lengths[,yPos := ceiling(cumsum(length)/length[1])]
yPos = 1
xStart = 0
y = 1
refLength = lengths$length[1] + 2000
cs = refLength +1
for(row in 2:nrow(lengths)) {
  l = lengths$length[row] + 2000
  if(cs + l > refLength) {
    cs = 0
    y = y + 1
  } 
  yPos = c(yPos, y)
  xStart = c(xStart, cs)
  cs = cs + l
}


contigCoords = data.table(lengths[,.(contig, length)], xStart, yPos)

genes[, c("yPos", "xStart") := contigCoords[match(Contig, contig), .(yPos, xStart)]]
cds[, c("yPos", "xStart") := contigCoords[match(scaffold, contig), .(yPos, xStart)]]
sdr[, c("yPos", "xStart") := contigCoords[match(CHROM, contig), .(yPos, xStart)]]
coding[, c("yPos", "xStart") := contigCoords[match(scaffold, contig), .(yPos, xStart)]]

margin = 0.2


pdf("../fig3.pdf", width = 9, height = 5)

with(contigCoords, plot(range(0, (length+xStart)/1000), range(margin, yPos + margin), ylab = "",
     xlab = "Position in contig (kbp)", bty = "n", type = "n", las = 1, yaxt = "n"))
with(contigCoords, rect(xStart/1000, yPos - margin, (length + xStart)/1000, yPos + margin, border = NA, col = "lightgrey"))

with(sdr[!is.na(yPos) & rec == 4L & !POS %in% coding$pos], segments((POS + xStart)/1000, yPos - margin, y1 = yPos + margin, col = "forestgreen", lend = 1))

with(coding[!is.na(yPos)], segments((pos+xStart)/1000, yPos - margin, y1 = yPos + margin, 
                                    col = ifelse(refAA == altAA, "royalblue", "red"), lend = 1))

with(genes[!is.na(yPos)], segments((Start + xStart)/1000, yPos + ifelse(Strand == "+", margin, -margin), (End + xStart)/1000, lend = 1, lwd = 2, col = "gray40"))
with(cds[!is.na(yPos)], segments((start + xStart)/1000, yPos + ifelse(strand == "+", margin, -margin), (end + xStart)/1000, lwd = 6, lend = 1))

with(genes[grep("Androgenic|Disintegrin", Product)], text((End+Start)/2000, yPos + ifelse(Strand == "+", margin, -margin), labels = c("agh","adam12"), cex = 0.8, pos = ifelse(Strand == "+", 3, 1)))

dev.off()
