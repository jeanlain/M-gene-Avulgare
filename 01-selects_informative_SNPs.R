##%######################################################%##
#                                                          #
####       Stage 1, processing parental genotypes       ####
####               for both families                     ###
#                                                          #
##%######################################################%##

library(data.table)

# Read the files using fread() and set column names
maf <- fread("parents_BF_first_family_angsd_dc.mafs", 
             col.names = c("chromo", "position", "major", "minor", "unknownEM", "pu_EM", "nInd"))
beagle <- fread("parents_BF_first_family_angsd_dc.beagle", 
                col.names = c("marker", "allele1", "allele2", "ind0_1", "ind0_2", "ind0_3", "ind1_1", "ind1_2", "ind1_3"))
counts <- fread("parents_BF_first_family_angsd_dc.counts")
pos <- fread("parents_BF_first_family_angsd_dc.pos")

# Combine the data tables
all <- cbind(beagle, maf, pos, counts) 

# Create a prototype VCF data table
proto_vcf <- all[, .(allele1, allele2, ind0_1, ind0_2, ind0_3, ind1_1, ind1_2, ind1_3, 
                     chromo, position, major, minor, unknownEM, pu_EM, nInd, chr, pos, 
                     totDepth, ind0_A, ind0_C, ind0_G, ind0_T, ind1_A, ind1_C, ind1_G, ind1_T)]

# Calculate total depths for individuals
proto_vcf[, totDepth_ind0 := ind0_A + ind0_C + ind0_G + ind0_T]
proto_vcf[, totDepth_ind1 := ind1_A + ind1_C + ind1_G + ind1_T]

# Filter data for male heterozygotes and female homozygotes
male_hetero <- proto_vcf[ind1_2 > 0.95]
male_hetero_female_homo <- male_hetero[ind0_1 > 0.95]
male_hetero_female_homo[, `:=`(pos = NULL, chr = NULL)]  # Remove 'pos' and 'chr' columns

# Increment allele counts for female homozygotes
male_hetero_female_homo[, `:=`(allele1 = allele1 + 1L, allele2 = allele2 + 1L)]

# Extract allele matrices for parents
mat_mother <- male_hetero_female_homo[, .(ind0_A, ind0_C, ind0_G, ind0_T)]
mat_father <- male_hetero_female_homo[, .(ind1_A, ind1_C, ind1_G, ind1_T)]

# Count alleles for each base
mother_major_count <- mat_mother[cbind(seq_len(nrow(mat_mother)), male_hetero_female_homo$allele1)]
mother_minor_count <- mat_mother[cbind(seq_len(nrow(mat_mother)), male_hetero_female_homo$allele2)]
father_major_count <- mat_father[cbind(seq_len(nrow(mat_father)), male_hetero_female_homo$allele1)]
father_minor_count <- mat_father[cbind(seq_len(nrow(mat_father)), male_hetero_female_homo$allele2)]

# Create a data table with allele counts
dt <- data.table(mother_major_count, mother_minor_count, father_major_count, father_minor_count)

# Merge data and compute totals
fusion_td_male_female <- cbind(male_hetero_female_homo, dt)
fusion_td_male_female[, `:=`(tot_maj_min_father = father_major_count + father_minor_count, 
                             tot_maj_min_mother = mother_major_count + mother_minor_count)]

# Filter data based on quantiles
depth_quantiles <- fusion_td_male_female[, quantile(tot_maj_min_mother, probs = 0.95, na.rm = TRUE)]
fusion_td_male_female <- fusion_td_male_female[tot_maj_min_mother < depth_quantiles & tot_maj_min_mother > 5L]

# Remove unnecessary columns
fusion_td_male_female[, `:=`(unknownEM = NULL, pu_EM = NULL, totDepth = NULL, 
                             ind0_A = NULL, ind0_C = NULL, ind0_G = NULL, ind0_T = NULL, 
                             ind1_A = NULL, ind1_C = NULL, ind1_G = NULL, ind1_T = NULL)]
fusion_td_male_female[, `:=`(ind0_2 = NULL, ind0_3 = NULL, ind1_1 = NULL, ind1_3 = NULL, 
                             major = NULL, minor = NULL, nInd = NULL, 
                             tot_maj_min_father = NULL, tot_maj_min_mother = NULL)]

# Reorganize and rename columns
fusion_td_male_female <- fusion_td_male_female[, .(allele1, allele2, ind0_1, ind1_2, chromo, position, 
                                                   totDepth_ind0, totDepth_ind1, mother_major_count, 
                                                   mother_minor_count, father_major_count, father_minor_count)]
setnames(fusion_td_male_female, c("allele1", "allele2", "ind0_1", "ind1_2", "CHROM", "POS", 
                                  "mother.DP", "father.DP", "A2_mother_count", "A1_mother_count", 
                                  "A2_father_count", "A1_father_count"))

# Write the results to a file
fwrite(fusion_td_male_female, "snp_informatifs_angsd.txt")
