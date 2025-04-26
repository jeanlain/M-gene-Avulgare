setwd("/Users/jean/Library/CloudStorage/OneDrive-UniversitédePoitiers/recherche/symchrosex/M_gene")
source("scripts_SNPs/WZ_functions.R")


beagle =fread("0-F0/fam1_2_mother_father_remapped_gatk.beagle", col.name = c("contig","allele1","allele2","fam_1_mo_1","fam_1_mo_2","fam_1_mo_3","fam_1_fa_1","fam_1_fa_2","fam_1_fa_3","fam_2_mo_1","fam_2_mo_2","fam_2_mo_3","fam_2_fa_1","fam_2_fa_2","fam_2_fa_3"))

counts  = fread("0-F0/fam1_2_mother_father_remapped_gatk.counts", col.names = c("fam_1_mo_A","fam_1_mo_C","fam_1_mo_G","fam_1_mo_T","fam_1_fa_A","fam_1_fa_C","fam_1_fa_G","fam_1_fa_T","fam_2_mo_A","fam_2_mo_C","fam_2_mo_G","fam_2_mo_T","fam_2_fa_A","fam_2_fa_C","fam_2_fa_G","fam_2_fa_T","va"))

pos= fread("0-F0/fam1_2_mother_father_remapped_gatk.pos")


proto_vcf = data.table(pos, beagle, counts)

proto_vcf = proto_vcf[,.(chr,pos,allele1,allele2,fam_1_fa_2,fam_2_fa_2,fam_1_fa_A,fam_1_fa_C,fam_1_fa_G,fam_1_fa_T,fam_2_fa_A,fam_2_fa_C,fam_2_fa_G,fam_2_fa_T)]

# computes sequencing depth on fathers
proto_vcf = proto_vcf[,totDepth_fam_1_fa := fam_1_fa_A+fam_1_fa_C+fam_1_fa_G+fam_1_fa_T]
proto_vcf = proto_vcf[,totDepth_fam_2_fa := fam_2_fa_A+fam_2_fa_C+fam_2_fa_G+fam_2_fa_T]

depthQuantiles = proto_vcf[, .(fam1 = quantile(totDepth_fam_1_fa, 0.95), fam2 = quantile(totDepth_fam_2_fa, 0.95))]

proto_vcf[, hetFa1 := fam_1_fa_2 >= 0.99 & totDepth_fam_1_fa < depthQuantiles$fam1 & totDepth_fam_1_fa >= 5L]
proto_vcf[, hetFa2 := fam_2_fa_2 >= 0.99 & totDepth_fam_2_fa < depthQuantiles$fam2 & totDepth_fam_2_fa >= 5L]

snpCount = proto_vcf[, .(oneFather = sum(hetFa1 | hetFa2), bothFathers = sum(hetFa1 & hetFa2)), by = chr]


writeT(snpCount, "tables/heteroSNPs_byContigP99q95depth5.txt")
writeT(depthQuantiles, "tables/depthQuantileFathers.txt")

