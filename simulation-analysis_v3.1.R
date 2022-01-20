# Author: Andrija Sente
# 
# Differential GABA-A receptor assembly diversifies structures and signaling
#
# Andrija Sente1*, Rooma Desai2, Katerina Naydenova1, Tomas Malinauskas3, 
# Youssef Jounaidi2, Jonas Miehling1, Xiaojuan Zhou2, Simonas Masiulis1,4, 
# Steven W. Hardwick5, Dimitri Y. Chirgadze5, Keith W. Miller2*, A. Radu Aricescu1*
#
# 1 - MRC Laboratory of Molecular Biology, Francis Crick Avenue, Cambridge, CB2 0QH, UK.
# 2 - Department of Anesthesia, Critical Care and Pain Medicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA, USA. 
# 3 - Division of Structural Biology, Wellcome Centre for Human Genetics, University of Oxford, Roosevelt Drive, Oxford, OX3 7BN, UK.
# 4 - Current address: Materials and Structural Analysis Division, Thermo Fisher Scientific, Achtseweg Noord, Eindhoven, 5651 GG, Netherlands.
# 5 - Department of Biochemistry, University of Cambridge, Tennis Court Road, Cambridge, CB2 1GA, UK.
#
# * - Correspondence to: asente@mrc-lmb.cam.ac.uk or k_miller@helix.mgh.harvard.edu or radu@mrc-lmb.cam.ac.uk

# The simulated data used in the paper is available at the following link:
# ftp://ftp.mrc-lmb.cam.ac.uk/pub/asente/differential-assembly/
  
# =============================================================================
# A note about nomenclature:
# Interface matrix:
#    a b  d
# a aa ab ad
# b ba bb bd
# d da db dd

# First subunit in the pairwise interface is principal, second is complementary
# ab = a+/b-
# All receptor arrangements are read COUNTERCLOCKWISE.

require(ggplot2)
library(pheatmap)
library(RColorBrewer)
require(reshape2)
require(dplyr)

setwd('~/your-directory/')

# Table of contents:

# 1. Preparing the data
#		1.1 Generating a receptor subtype lookup table
#		1.2 Reading the simulation data
# 2. Searching for conditions compatible with alpha4/beta3/delta cryoEM data
# 		2.1 Find the requried conditions
#		2.2 Plots
# 3. Searching for conditions compatible with alpha4/beta3/gamma2 cryoEM data
# 		3.1 Find the requried conditions
#		3.2 Plots
# 4. Searching for mutually consistent conditions identified in (2) and (3)
#		This data is used in Figure 4 and Supplementary Figure 4
#		4.1 Generating a random subset for data in Fig.4 and SI Fig. 4
#		4.2 Plots with all data
#		4.3 Plots with subsets of data
#			4.3.1 - Figure 4b & Supplementary Figure 4a
#			4.3.2 - Supplementary Figure 4b
#			4.3.3 - Supplementary Figure 4c
#			4.4.4 - Saving panels for Fig. 4b & Supplementary Figure 4a-c
#			4.3.5 - Supplementary Figure 4d
#			4.3.6 - Supplementary Figure 4e
#			4.3.7 - Saving panels for Supplementary Figure 4d-e
# 5. Can we alter receptor subtype expression by alterning only the subunit abundance?
#		5.1 a4b3d case:Supplementary Figure 5b - all panels
#		5.2 a4b3d case:Saving panels for Supplementary Figure 5b
#		5.3 a4b3g2 case:Supplementary Figure 5c - all panels
#		5.4 a4b3g2 case:Saving panels for Supplementary Figure 5c
# 6. Data for Supplementary Figure 6a 
# 7. Identifying conditions compatible with a1b3g2 cryo-EM data
# 		7.1 Data for Supplementary Figure 5d
# 		7.2 Figure panels for Supplementary Figure 5d 
# 		7.3 Saving figure panels for Supplementary Figure 5d 
# 		7.4 Sample size for Supplementary Figure 5:
# 8. Boxplots and violin plots for Supplementary Fig. 5b-d
# 		8.1 Preparing data
# 		8.2 Supplementary Figure 5b - violin plot for interface likelihoods:
# 		8.3 Supplementary Figure 5d - violin plot for interface likelihoods:
# 		8.4 Supplementary Figure 5b,d - preparing data for abundance boxplots:
# 		8.5 Supplementary Figure 5b - abundance boxplots (a4b3d):
# 		8.6 Mean abundances identified
# 		8.7 Saving boxplots and violin plots for Supplementary Fig. 5b-d


# ==============================================================================

# Some common useful  receptors:
# ababb= Num 20 grep(pattern = 'ababb', x = lookup$combined)
# dbaba= Num 22 (this is the proper ABG 2:2:1) 
#         grep(pattern = 'dbaba', x = lookup$combined)
# dbabb= Num. 29 grep(pattern = 'dbabb', x = lookup$combined)
# adbbb = Num 27 (this is the ABBBD type) 
#         grep(pattern = 'adbbb', x = lookup$combined)
# dbdbb= Num. 47 grep(pattern = 'dbdbb', x = lookup$combined)
# bdbbb = Num. 45 grep(pattern = 'bdbbb', x = lookup$combined)
# b-homomer = Num. 44 grep(pattern = 'bbbbb', x = lookup$combined)



# ================================= PART 1 =====================================

# 1.1 Generating a lookup table for the 51 receptor subtypes that can assemble 
# 	from three subunits.

lookup <- read.table('receptor-numbering.txt', sep = ' ')
lookup$receptor <- apply( lookup[ , 1:5 ] , 1 , paste , collapse = "" )
lookup$V6 <- lookup$V6 + 1
lookup <- subset(lookup, select = c('receptor', 'V6'))
colnames(lookup) <- c('receptorType', 'receptorNumber')

# any  circular permutation is contained within the concatenated string - this
# enables easier searching later on
lookup$combined <- paste(lookup$receptorType, lookup$receptorType, sep = "")


# 1.2 Loading and annotating the simulated data
data <- read.table('simulation_aa0_all.txt', sep = ' ')
colnames(data) <- c('alpha', 'beta', 'delta', 'aa', 'ab', 'ad', 
                    'ba', 'bb', 'bd', 'da', 'db', 'dd', lookup$receptorType)

# removing the first 12 columns of the simulation data
# these are initiated interface abundances (1-3) and affinities (4-12)

just_data <- data[,13:ncol(data)]
interfaces <- log10(data[,4:12])
interfaces[interfaces == '-Inf'] <- -4
interfaces[interfaces < -1 ] <- -1

abundance <- data[,1:3]




# ================================= PART 2 =====================================

# Searching for conditions that have ADBBB (counterclockwise) and BBBBD subtypes
# predominant from the population of receptors that contain the delta subunit, 
# but don't have other "solvable" arrangements

# 2.1 Find the requried conditions

# ----- (1) Find all solvable receptors:
# We can solve any receptor that has (A) a delta-subunit and (B) a b/b interface

INDEX_delta_containing <- grep(pattern = 'd', x = lookup$combined)
INDEX_BBinterface <- grep(pattern = 'bb', x = lookup$combined)

# To identify SOLVABLE but NOT OCCURING receptors, find intersection between the 
# 	two variables above
# Remove the solved b3d receptor (BBBBD = 45) and a4b3d (BBBAD = 27)

INDEX_notOccurring_a4b3d <- intersect(INDEX_delta_containing, INDEX_BBinterface)

INDEX_notOccurring_a4b3d <- 
  INDEX_notOccurring_a4b3d[ ! INDEX_notOccurring_a4b3d %in% c('27', '45')]

# ----- (2) Find sums of all receptors with:
# (A) - delta subunits
# (B) - receptors we observe (abbbd and bbbbd)
sumAllDelta_a4b3d <- apply(just_data[,INDEX_delta_containing], 1, sum) #purified
sumOurReceptors_a4b3d <- apply(just_data[,c(27, 45)], 1, sum) # solved receptors

# ----- (3) Find rows in the data that satisfy the conditions:
# (A) - Obsesrved receptors are 50% or more of purified receptors
# (B) - Delta-containing receptors take up >50% of all expressed receptors ***
# (C) - Observed receptors are above noise level (>33)
# (D) - Not-observed receptors are below the noise level (<33)
# (E) - Both BBBAD and BBBBD are abundant (>1/3 of purified receptors each) ***

rowOfConditions_a4b3d <- which(   
  
  # (A) - Obsesrved receptors are 50% or more of purified receptors
  (sumOurReceptors_a4b3d > sumAllDelta_a4b3d/2) & 
    
    # (B) - Delta-containing receptors take up >50% of all expressed receptors ***
    (sumAllDelta_a4b3d > 300)  & 
    
    # (C) - Observed receptors are above noise level (>33)
    (just_data[,27] > 33 ) & # BBBAD
    (just_data[,45] > 33 ) & # BBBBD
    
    # (D) - Not-observed receptors are below the noise level (<33)
    # (just_data[,INDEX_notOccurring_a4b3d] < 33 ) &
    (just_data[,12] < 33 ) &
    (just_data[,15] < 33 ) &
    (just_data[,25] < 33 ) &
    (just_data[,26] < 33 ) &
    (just_data[,28] < 33 ) &
    (just_data[,29] < 33 ) &
    (just_data[,36] < 33 ) &
    (just_data[,38] < 33 ) &
    (just_data[,40] < 33 ) &
    (just_data[,46] < 33 ) &
    (just_data[,47] < 33 ) &
    (just_data[,48] < 33 ) &
    
    
    # (E) - Both BBBAD and BBBBD are abundant (>1/3 of purified receptors each)**
    (just_data[,27] > sumOurReceptors_a4b3d/3) & # BBBAD
    (just_data[,45] > sumOurReceptors_a4b3d/3) ) # BBBBD

# length(rowOfConditions_a4b3d) = 211100



# 2.2 Plots:
# NOTE: these plots are massive

heatmapIdentifiedConditions_a4b3d <- pheatmap(
  mat = just_data[rowOfConditions_a4b3d, ], 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')

heatmapConditionsAbundance_a4b3d <- pheatmap(
  mat = abundance[rowOfConditions_a4b3d, ], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')

heatmapConditionsInterfaces_a4b3d <- pheatmap(
  mat = interfaces[rowOfConditions_a4b3d, ], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')




# ================================= PART 3 =====================================

# 3. Searching for conditions compatible with alpha4/beta3/gamma2 cryoEM data

# Searching for conditions in which the predominant receptor subtypes, among the
# purifiable receptors, are the BBGBA (counterclockwise) and BBGBG subtypes, but 
# other "solvable" arrangements are kept at minimum.

# 3.1 Find the requried conditions

# ----- (1) Find all solvable receptors:
# We can solve any receptor that has (A) a delta-subunit and (B) a b/b interface
INDEX_gamma_containing <- grep(pattern = 'd', x = lookup$combined)
INDEX_BBinterface <- grep(pattern = 'bb', x = lookup$combined)

# To identify SOLVABLE but NOT OCCURING receptors, find intersection between the
# two variables above
# Remove the solved b3g receptor (BBGBG = 47) and a4b3g (BBGBA = 29)
INDEX_notOccurring_gamma <- intersect(INDEX_gamma_containing, INDEX_BBinterface)

INDEX_notOccurring_gamma <- 
  INDEX_notOccurring_gamma[ ! INDEX_notOccurring_gamma %in% c('29', '47')]

# ----- (2) Find sums of all receptors with:
# (A) - gamma subunits
# (B) - receptors we observe (BBGBG and BBGBA)
sumAllGamma <- apply(just_data[,INDEX_gamma_containing], 1, sum) # purified 
sumOurReceptors_gamma <- apply(just_data[,c(29, 47)], 1, sum) # solved 

# ----- (3) Find rows in the data that satisfy the conditions:
# (A) - Observed receptors are 50% or more of purified receptors
# (B) - Gamma-containing receptors take up >50% of all expressed receptors ***
# (C) - Observed receptors are above noise level (>33)
# (D) - Not-observed receptors are below the noise level (<33)
# (E) - Both BBGBA and BBGBG are abundant (>1/3 of purified receptors each) ***

rowOfConditions_gamma <- which( 
  
  # (A) - Obsesrved receptors are 50% or more of purified receptors
  (sumOurReceptors_gamma > sumAllGamma/2) & 
    
    # (B) - Delta-containing receptors take up >50% of all expressed receptors ***
    (sumAllGamma > 300)  & 
    
    # (C) - Observed receptors are above noise level (>33)
    (just_data[,29] > 33 ) & # BBGBA
    (just_data[,47] > 33 ) & # BBGBG
    
    # (D) - Not-observed receptors are below the noise level (<33)
    # (just_data[,INDEX_notOccurring_gamma] < 33) &
    
    (just_data[,12] < 33 ) &
    (just_data[,15] < 33 ) &
    (just_data[,25] < 33 ) &
    (just_data[,26] < 33 ) &
    (just_data[,27] < 33 ) &
    (just_data[,28] < 33 ) &
    (just_data[,36] < 33 ) &
    (just_data[,38] < 33 ) &
    (just_data[,40] < 33 ) &
    (just_data[,45] < 33 ) &
    (just_data[,46] < 33 ) &
    (just_data[,48] < 33 ) &
    
    # (E) - Both BBGBA and BBGBG are abundant (>1/3 of purified receptors each)**
    (just_data[,29] > sumOurReceptors_gamma/3) & # BBGBA
    (just_data[,47] > sumOurReceptors_gamma/3) ) # BBGBG

# length(rowOfConditions_gamma) = 1524

# 3.2 Plots:
# NOTE: These plots are massive

heatmapIdentifiedConditions_gamma <- pheatmap(
  mat = just_data[rowOfConditions_gamma, ], 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')

heatmapConditionsAbundance_gamma <- pheatmap(
  mat = abundance[rowOfConditions_gamma, ], 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')

heatmapConditionsInterfaces_gamma <- pheatmap(
  mat = interfaces[rowOfConditions_gamma, ], 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')


# ================================= PART 4 =====================================

# In this part, we search for cases among the alpha4/beta3/delta (section 2) and
# alpha4/beta3/gamma2 (section 3) identified conditions in which AA, AB, BA and BB 
# interfaces have the same values

# ----- (1) Take the data that satisfies criteria given above
data_gamma <- data[rowOfConditions_gamma, ] # data satisfying a4b3g2 cryoEM
data_gamma$dataset <- 'gamma'
data_delta <- data[rowOfConditions_a4b3d, ] # data satisfying a4b3d cryoEM
data_delta$dataset <- 'delta'

# ----- (2) Concatenating interface columns to enable easier searching
data_gamma$first_four <- apply( 
  data_gamma[ ,c('aa', 'ab', 'ba', 'bb') ] , 1 , paste , collapse = "" )
data_delta$first_four <- apply( 
  data_delta[ ,c('aa', 'ab', 'ba', 'bb') ] , 1 , paste , collapse = "" )

# ----- (3) Search for matching parameters in "delta" and "gamma" conditions 
data_gamma <- subset(data_gamma, 
                     data_gamma$first_four %in% data_delta$first_four)
data_delta <- subset(data_delta, 
                     data_delta$first_four %in% data_gamma$first_four)

# ----- (4) Combine data and arrange based on dataset (delta/gamma)
data_combined <- rbind(data_gamma, data_delta)
data_combined <- arrange(data_combined, dataset)

data_combined$dataset <- NULL
data_combined$first_four <- NULL

# ----- (5 - FINAL) Combined delta/gamma INTERFACE data
interfacesCombined <- log10(data_combined[,4:12])
interfacesCombined[interfacesCombined == '-Inf'] <- -4
interfacesCombined[interfacesCombined < -1 ] <- -1

# ----- (6 - FINAL) Combined delta/gamma ABUNDANCE data
abundance <- data_combined[,1:3]

# ----- (7 - FINAL) Combined delta/gamma SUBTYPE EXPRESSION data
just_data_combined <- data_combined[,13:ncol(data_combined)]

# nrow(data_gamma) = 72
# nrow(data_delta) = 2967
# These are the conditions shown in Figure 4 and Supplementary Figure 4.


# 4.1
# Take 20 random rows from mutually consistent a4b3d and a4b3g2 data compatible 
# with cryoEM (data_combined) and make panels for Figure 4 and SI Figure 4

data_gamma_sample <- data_gamma[sample(nrow(data_gamma), 20), ]
data_delta_sample <- data_delta[sample(nrow(data_delta), 20), ]

# ----- (4) Combine data and arrange based on dataset (delta/gamma)
data_combined_sample <- rbind(data_gamma_sample, data_delta_sample)
data_combined_sample <- arrange(data_combined_sample, dataset)

data_combined_sample$dataset <- NULL
data_combined_sample$first_four <- NULL

# ----- (5 - FINAL) Combined delta/gamma INTERFACE data
interfacesCombined_sample <- log10(data_combined_sample[,4:12])
interfacesCombined_sample[interfacesCombined_sample == '-Inf'] <- -4
interfacesCombined_sample[interfacesCombined_sample < -1 ] <- -1

# ----- (6 - FINAL) Combined delta/gamma ABUNDANCE data
abundance_sample <- data_combined_sample[,1:3]

# ----- (7 - FINAL) Combined delta/gamma ABUNDANCE data
just_data_combined_sample <- data_combined_sample[,13:ncol(data_combined_sample)]
just_data_combined_sample <- 
  just_data_combined_sample[, -which(colSums(just_data_combined_sample)==0)]




# 4.2
# Plotting interface affinity, subunit abundance and receptor subtype expression
# for all conditions compatible with cryoEM data.

combined_data <- pheatmap(
  mat = just_data_combined, 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')

combined_interfaces <- pheatmap(
  mat = interfacesCombined, 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')

combined_abundances <- pheatmap(
  mat = abundance, 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45')



# 4.3
# These plots are made from 20 randomly selected conditions from 4.2
# This generates panels for Figure 4 and Supplementary Figure 4.

# 4.3.1 - Figure 4b & Supplementary Figure 4a
combined_data <- pheatmap(
  mat = just_data_combined_sample, 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45',
  show_rownames = F)

# 4.3.2 - Supplementary Figure 4b
combined_interfaces <- pheatmap(
  mat = interfacesCombined_sample, 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100),
  angle_col = '45',
  show_rownames = F)

# 4.3.3 - Supplementary Figure 4c
combined_abundances <- pheatmap(
  mat = abundance_sample, 
  cluster_cols = F, cluster_rows = F, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
  angle_col = '45',
  show_rownames = F)



# 4.4.4 - Saving figure panels

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/Subtypes20.pdf', 
#     width = 4, height = 2, onefile = F) 
# print(combined_data)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/Interface20.pdf', 
#     width = 1.6, height = 2, onefile = F) 
# print(combined_interfaces)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/Abundance20.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(combined_abundances)
# dev.off()



# 4.4.5 - Supplementary Figure 4d
# Generating a violin plot of ddG of subunit binding (kcal/mol)

require(tidyr)

# Convert the interface likelihood to actual ddG
# kT=RT/N(Avogadro) at 310 K (37 deg C)
# To convert the likelihood into relative energy, use the formula:
# p ~ e^(-epsilon/(kT))

interfacesGamma_boxplot <- -log(data_gamma[,4:12])*0.593
interfacesDelta_boxplot <- -log(data_delta[,4:12])*0.593

# prepare the data: 
interfacesGamma_boxplot <- pivot_longer(
  data = interfacesGamma_boxplot,
  cols = colnames(interfacesGamma_boxplot)[1:9],
  names_to = "interface",
)

interfacesGamma_boxplot$receptor <- 'gamma'

interfacesDelta_boxplot <- pivot_longer(
  data = interfacesDelta_boxplot,
  cols = colnames(interfacesDelta_boxplot)[1:9],
  names_to = "interface",
)

interfacesDelta_boxplot$receptor <- 'delta'

interfacesCombined_boxplot <- rbind(interfacesGamma_boxplot,
                                    interfacesDelta_boxplot)

interfacesCombined_boxplot <- subset(interfacesCombined_boxplot, 
                                     interfacesCombined_boxplot$interface != 'aa')

# plot the violing plot
interface_boxplot <- ggplot(data = interfacesCombined_boxplot, 
                            aes(x = interfacesCombined_boxplot$interface, 
                                y = interfacesCombined_boxplot$value, 
                                fill = interfacesCombined_boxplot$receptor,
                                color = interfacesCombined_boxplot$receptor)) +
  geom_violin(scale='area') +
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey", size=0.25)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_line(color='grey', size=0.25),
        legend.position = "none",
        text = element_text(size=7)) #+
# scale_fill_manual(values=c("#a6611a", "#018571"))+
# scale_color_manual(values=c("#a6611a", "#018571"))

# scale_fill_manual(values=c("#d7191c", "#2c7bb6")) +
# scale_color_manual(values=c("#d7191c", "#2c7bb6")) 
# scale_fill_manual(values=c("#0571b0", "#ca0020")) +
# scale_color_manual(values=c("#0571b0", "#ca0020")) 

interface_boxplot 

# interface_boxplot <- ggplot(data = interfacesCombined_boxplot, 
#                             aes(x = interfacesCombined_boxplot$interface, 
#                                 y = interfacesCombined_boxplot$value, 
#                                 alpha = interfacesCombined_boxplot$interface)) +
#   geom_boxplot(notch=TRUE)
# interface_boxplot



# 4.3.6 - Supplementary Figure 4e
# Creates nothced box plots for comparison of subunit abundances

abundancesDelta_boxplot <- data_delta[,1:3]
abundancesGamma_boxplot <- data_gamma[,1:3]

abundancesDelta_boxplot <- pivot_longer(
  data = abundancesDelta_boxplot,
  cols = colnames(abundancesDelta_boxplot)[1:3],
  names_to = "interface",
)

abundancesGamma_boxplot <- pivot_longer(
  data = abundancesGamma_boxplot,
  cols = colnames(abundancesGamma_boxplot)[1:3],
  names_to = "interface",
)

abundancesDelta_boxplot$receptor <- 'delta'
abundancesGamma_boxplot$receptor <- 'gamma'

abundance_boxplot <- rbind(abundancesGamma_boxplot,abundancesDelta_boxplot)

# create the boxplot:
abundanceBoxplot <- ggplot(data = abundance_boxplot, 
                           aes(x = abundance_boxplot$interface, 
                               y = abundance_boxplot$value,
                               fill = abundance_boxplot$receptor,
                               color = abundance_boxplot$receptor)) +
  # ),color='gray') +
  # geom_violin(scale='area') +
  # geom_boxplot(notch=TRUE, outlier.size=0.2, lwd=0.25)+
  geom_boxplot(notch=TRUE, outlier.shape = NA) +
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey", size=0.25)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_line(color='grey', size=0.25),
        legend.position = "none",
        text = element_text(size=7)) +
  theme(legend.position = "none") #+
# scale_fill_manual(values=c("#0571b0", "#ca0020")) +
# scale_color_manual(values=c("#0571b0", "#ca0020")) 

#+
# scale_fill_manual(values=c("#fc9272", "#de2d26"))
abundanceBoxplot

# Means for a4b3d
# mean(subset(abundancesDelta_boxplot, interface == 'alpha')$value) = 0.02804335
# mean(subset(abundancesDelta_boxplot, interface == 'beta')$value) = 0.689431
# mean(subset(abundancesDelta_boxplot, interface == 'delta')$value) = 0.2825257
# 
# Means for a4b3g2
# mean(subset(abundancesGamma_boxplot, interface == 'alpha')$value) = 0.05090304
# mean(subset(abundancesGamma_boxplot, interface == 'beta')$value) = 0.02161769
# mean(subset(abundancesGamma_boxplot, interface == 'delta')$value) = 0.9274793


# 4.3.7 - save the figure panels for Supplementary Figure 4d-e

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/InterfaceBoxplot.pdf', 
#     width = 1.3, height = 2, onefile = F) 
# print(interface_boxplot)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/AbundanceBoxplot.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(abundanceBoxplot)
# dev.off()










# ================================= PART 5 =====================================
# 
# Here, we ask a question if it is possible to obtain a receptor with the 
# canonical babad (label #22) arrangement simply by altering subunit abundances.
# To investigate this, we keep the interface affinities the same as those found
# in the a4b3d (bbabd) case and search for the conditions in which babad is the
# most prevalent.
# The only other condition is that, among purifiable receptors, dbaba is the 
# most prevalent. Other "solvable" arrangements are kept at minimum.

# 5.1 - Generates data and panels for Supplementary Figure 5b
# This is analysis for alpha4/beta3/delta

data_babad <- subset(data, data$aa %in% data_delta$aa &
                       data$ab %in% data_delta$ab &
                       data$ad %in% data_delta$ad &
                       data$ba %in% data_delta$ba &
                       data$bb %in% data_delta$bb &
                       data$bd %in% data_delta$bd &
                       data$da %in% data_delta$da &
                       data$db %in% data_delta$db &
                       data$dd %in% data_delta$dd)


require(plyr)
# these two lines of code return the matched rows from the first dataset
data_babad <- plyr::match_df(data_babad, data_delta, 
                             c("aa", "ab", "ad", "ba", "bb", "bd", "da", "db", "dd"))
data_babad <- plyr::match_df(data_babad, data_gamma, 
                             c("aa", "ab", "ba", "bb"))

# subsetting is correct:
# > nrow(unique(data_babad[,4:12]))
# [1] 505
# > nrow(unique(data_delta[,4:12]))
# [1] 505


# Creating useful dataframes:

just_data_babad <- data_babad[,13:ncol(data_babad)]
interfaces_babad <- log10(data_babad[,4:12])
interfaces_babad[interfaces_babad == '-Inf'] <- -4
interfaces_babad[interfaces_babad < -1 ] <- -1
abundance_babad <- data_babad[,1:3]

# ----- (1) Find purifiable receptors:
# We can purify any receptor that has a delta subunit
INDEX_delta_babad <- grep(pattern = 'd', x = lookup$combined) 

# ----- (2) Find sums of all receptors with:
# (A) - delta_babad subunits
# (B) - receptors we observe (BBGBG and BBGBA)
sumAlldelta_babad <- apply(just_data_babad[,INDEX_delta_babad], 1, sum) #purified 

# ----- (3) Find rows in the data that satisfy the conditions:
# (A) - Observed receptors are 50% or more of purified receptors
# (B) - delta-containing receptors take up >50% of all expressed receptors ***

rowOfConditions_babad <- which( 
  
  # (A) - Obsesrved receptors are 50% or more of purified receptors
  (data_babad$dbaba > sumAlldelta_babad/2) & 
    
    # (B) - our receptor of interest is >60% of all expressed receptors ***
    (just_data_babad[,c(22)] > 600)  #& 
  
)

# length(rowOfConditions_babad) = 4853


# Panels for Supplementary Figure 5b:

babad_abundance <- pheatmap(
  mat = abundance_babad[rowOfConditions_babad,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

babad_interface <- pheatmap(
  mat = interfaces_babad[rowOfConditions_babad,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

babad_expression <- pheatmap(
  mat = just_data_babad[rowOfConditions_babad,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

# 5.2 Saving panels for Supplementary Figure 5b

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babad_expression.pdf', 
#     width = 4, height = 2, onefile = F) 
# print(babad_expression)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babad_interface.pdf', 
#     width = 1.6, height = 2, onefile = F) 
# print(babad_interface)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babad_abundance.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(babad_abundance)
# dev.off()






# 5.3 - Generates data and panels for Supplementary Figure 5c
# This is analysis for alpha4/beta3/gamma2
# The interfaces are kept the same as in data_gamma

data_babag <- plyr::match_df(data, data_gamma, 
                             c("aa", "ab", "ad", "ba", "bb", "bd", "da", "db", "dd"))
data_babag <- plyr::match_df(data_babag, data_delta, 
                             c("aa", "ab", "ba", "bb"))

# subsetting is correct:
# > nrow(unique(data_babag[,4:12]))
# [1] 22
# > nrow(unique(data_gamma[,4:12]))
# [1] 22


# Generate useful dataframes:

just_data_babag <- data_babag[,13:ncol(data_babag)]
interfaces_babag <- log10(data_babag[,4:12])
interfaces_babag[interfaces_babag == '-Inf'] <- -4
interfaces_babag[interfaces_babag < -1 ] <- -1
abundance_babag <- data_babag[,1:3]

# ----- (1) Find purifiable receptors:
# We can purify any receptor that has a gamma subunit
INDEX_gamma <- grep(pattern = 'd', x = lookup$combined)

# ----- (2) Find sums of all receptors with:
# (A) - gamma subunits
# (B) - receptors we observe (BBGBG and BBGBA)
sumAllGamma_babag <- apply(just_data_babag[,INDEX_gamma], 1, sum) # purified 

# ----- (3) Find rows in the data that satisfy the conditions:
# (A) - Observed receptors are 50% or more of purified receptors
# (B) - Gamma-containing receptors take up >60% of all expressed receptors ***

rowOfConditions_babag <- which( 
  
  # (A) - Obsesrved receptors are 50% or more of purified receptors
  (data_babag$dbaba > sumAllGamma_babag/2) & 
    
    # (B) - Delta-containing receptors take up >60% of all expressed receptors ***
    (just_data_babag[,c(22)] > 600) )

# length(rowOfConditions_babag)=3


# These are panels for Supplementary Figure 5d

babag_abundance <- pheatmap(
  mat = abundance_babag[rowOfConditions_babag,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0,
  show_colnames = F,
  legend = FALSE)

babag_interface <- pheatmap(
  mat = interfaces_babag[rowOfConditions_babag,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0,
  show_colnames = F,
  legend = FALSE)

babag_expression <- pheatmap(
  mat = just_data_babag[rowOfConditions_babag,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0,
  show_colnames = F,
  legend = FALSE)

# 5.4 Saving panels for Supplementary Figure 5d

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babag_expression.pdf', 
#     width = 4, height = 0.2, onefile = F) 
# print(babag_expression)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babag_interface.pdf', 
#     width = 1.6, height = 0.2, onefile = F) 
# print(babag_interface)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babag_abundance.pdf', 
#     width = 1, height = 0.2, onefile = F) 
# print(babag_abundance)
# dev.off()






# 6. Data for Supplementary Figure 6a 
# takes a random subset of 100k of simulated conditions and plots a heatmap. 

# expression_profiles_sample <- just_data[sample(nrow(just_data), 100), ]
# expression_profiles_sample_10k <- just_data[sample(nrow(just_data), 10000), ]
# expression_profiles_examples <- pheatmap(
#   mat = expression_profiles_sample_10k[,20:ncol(expression_profiles_sample_10k)], 
#   cluster_cols = F, cluster_rows = T, fontsize = 7, 
#   color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
#   angle_col = '45',
#   show_rownames = F,
#   treeheight_row = 0, 
#   treeheight_col = 0)
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/expression_sample10k.pdf', 
#     width = 5, height = 3, onefile = F) 
# print(expression_profiles_examples)
# dev.off()

expression_profiles_sample_100k <- just_data[sample(nrow(just_data), 100000), ]

expression_profiles__100k <- pheatmap(
  mat = expression_profiles_sample_100k[,20:ncol(expression_profiles_sample_100k)], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/expression_profiles__100k.pdf', 
#     width = 5, height = 3, onefile = F) 
# print(expression_profiles_examples)
# dev.off()








# ================================= PART 7 =====================================

# In this part, we search for conditions in which a1b3g2 (BABAG arrangement) 
# is the prevalent receptor, with the restraint that relevant interfaces are the
# same as in the a4b3g2 receptor.
# Namely, these are the b3/b3, b3/g2, g2/b3 and g2/g2 interfaces.
# Also, we impose restraints that this arrangement is the most prevalent and
# other solvable and purifiable receptors don't occur (as shown by cryo=EM)
# This is presented in the paper as Supplementary Figure 5d

# 7.1 Data for Supplementary Figure 5d 

# retain only those simulated conditions in which the relevant interface 
# affinities are the same as those in a4b3g2 simulations.

data_a1bg <- plyr::match_df(data, data_gamma, 
                            c("bb", "bd", "db", "dd"))

# Some useful dataframes:

just_data_a1bg <- data_a1bg[,13:ncol(data_a1bg)]
interfaces_a1bg <- log10(data_a1bg[,4:12])
interfaces_a1bg[interfaces_a1bg == '-Inf'] <- -4
interfaces_a1bg[interfaces_a1bg < -1 ] <- -1
abundance_a1bg <- data_a1bg[,1:3]

# ----- (1) Find purifiable receptors:
# We can purify any receptor that has a gamma subunit
INDEX_gamma <- grep(pattern = 'd', x = lookup$combined)
INDEX_ABinterface <- grep(pattern = 'ab', x = lookup$combined)
INDEX_notOccurring_a1bg <- intersect(INDEX_gamma, INDEX_ABinterface)
INDEX_notOccurring_a1bg <-
  INDEX_notOccurring_a1bg [ ! INDEX_notOccurring_a1bg %in% c('22')]

# ----- (2) Find sums of all receptors with:
# (A) - gamma subunits
# (B) - receptors we observe (BBGBG and BBGBA)
sumAllGamma_a1bg <- apply(just_data_a1bg[,INDEX_gamma], 1, sum) # purified 

# ----- (3) Find rows in the data that satisfy the conditions:
# (A) - Observed receptors are 50% or more of purified receptors
# (B) - Gamma-containing receptors take up >50% of all expressed receptors ***
# (C) - Observed receptors are above noise level (>33)
# (D) - Not-observed receptors are below the noise level (<33)
# (E) - Both BBGBA and BBGBG are abundant (>1/3 of purified receptors each) ***

rowOfConditions_a1bg <- which( 
  
  # (A) - Obsesrved receptors are 50% or more of purified receptors
  (data_a1bg$dbaba > sumAllGamma_a1bg/2) & 
    
    # (B) - Delta-containing receptors take up >60% of all expressed receptors ***
    (just_data_a1bg[,22] > 600)  & 
    
    # # (C) - Observed receptors are above noise level (>33)
    (just_data_a1bg[,22] > 33 )  &
    
    # # (D) - Not-observed receptors are below the noise level (<33)
    # # (just_data[,INDEX_notOccurring_gamma] < 33) &
    
    (just_data_a1bg[,5] < 33 ) &
    (just_data_a1bg[,10] < 33 ) &
    (just_data_a1bg[,11] < 33 ) &
    (just_data_a1bg[,14] < 33 ) &
    (just_data_a1bg[,15] < 33 ) &
    (just_data_a1bg[,16] < 33 ) &
    (just_data_a1bg[,21] < 33 ) &
    (just_data_a1bg[,23] < 33 ) &
    (just_data_a1bg[,25] < 33 ) &
    (just_data_a1bg[,26] < 33 ) &
    (just_data_a1bg[,29] < 33 ) &
    (just_data_a1bg[,30] < 33 ) &
    (just_data_a1bg[,34] < 33 ) &
    (just_data_a1bg[,36] < 33 ) &
    (just_data_a1bg[,37] < 33 ) &
    (just_data_a1bg[,40] < 33 ) &
    (just_data_a1bg[,41] < 33 )
)

# length(rowOfConditions_a1bg)=27581

# babad_combined <- cbind(abundance_babad[rowOfConditions_babad,],
#                         interfaces_babad[rowOfConditions_babad,],
#                         just_data_babad[rowOfConditions_babad,])
# babad_combined <- dplyr::arrange(babad_combined, alpha, beta, delta)


# 7.2 Figure panels for Supplementary Figure 5d 

a1bg_abundance <- pheatmap(
  mat = abundance_a1bg[rowOfConditions_a1bg,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

a1bg_interface <- pheatmap(
  mat = interfaces_a1bg[rowOfConditions_a1bg,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

a1bg_expression <- pheatmap(
  mat = just_data_a1bg[rowOfConditions_a1bg,], 
  cluster_cols = F, cluster_rows = T, fontsize = 7, 
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  angle_col = '45',
  show_rownames = F,
  treeheight_row = 0, 
  treeheight_col = 0)

# 7.3 Saving figure panels for Supplementary Figure 5d 

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/a1bg_expression.pdf', 
#     width = 6, height = 4, onefile = F) 
# print(a1bg_expression)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/a1bg_interface.pdf', 
#     width = 3.2, height = 4, onefile = F) 
# print(a1bg_interface)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/a1bg_abundance.pdf', 
#     width = 2, height = 4, onefile = F) 
# print(a1bg_abundance)
# dev.off()

# 7.4 Sample size for Supplementary Figure 5:
# 5b - length(rowOfConditions_babad)=4853
# 5c - length(rowOfConditions_babag)=3
# 5d - length(rowOfConditions_a1bg)=27581




# ================================= PART 8 =====================================
# 8. Boxplots and violin plots for Supplementary Fig. 5b-d

# 8.1 Preparing data

require(tidyr)

interfaces_babad_boxplot <- -log(data_babad[rowOfConditions_babad,4:12])*0.593
interfaces_babag_boxplot <- -log(data_babag[rowOfConditions_babag,4:12])*0.593
interfaces_a1bg_boxplot <- -log(data_a1bg[rowOfConditions_a1bg,4:12])*0.593

# rotate the dataframes:

interfaces_babag_boxplot <- pivot_longer(
  data = interfaces_babag_boxplot,
  cols = colnames(interfaces_babag_boxplot)[1:9],
  names_to = "interface",
)

interfaces_babad_boxplot <- pivot_longer(
  data = interfaces_babad_boxplot,
  cols = colnames(interfaces_babad_boxplot)[1:9],
  names_to = "interface",
)

interfaces_a1bg_boxplot <- pivot_longer(
  data = interfaces_a1bg_boxplot,
  cols = colnames(interfaces_a1bg_boxplot)[1:9],
  names_to = "interface",
)


# 8.2 Supplementary Figure 5b - violin plot for interface likelihoods:

interfaces_babad_boxplot <- subset(interfaces_babad_boxplot,
                                   interfaces_babad_boxplot$interface != 'aa')
interface_babad_boxplot <- ggplot(data = interfaces_babad_boxplot, 
                                  aes(x = interfaces_babad_boxplot$interface, 
                                      y = interfaces_babad_boxplot$value, 
                                      color = "#9b87c4",
                                      fill = '#9b87c4')) +
  geom_violin(scale='area') +
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey", size=0.25)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_line(color='grey', size=0.25),
        legend.position = "none",
        text = element_text(size=7)) +
  scale_fill_manual(values=c("#9b87c4")) +
  scale_color_manual(values=c("#9b87c4"))

interface_babad_boxplot 


# 8.3 Supplementary Figure 5d - violin plot for interface likelihoods:

interfaces_a1bg_boxplot <- subset(interfaces_a1bg_boxplot,
                                  interfaces_a1bg_boxplot$interface != 'aa')
interface_a1bg_boxplot <- ggplot(data = interfaces_a1bg_boxplot, 
                                 aes(x = interfaces_a1bg_boxplot$interface, 
                                     y = interfaces_a1bg_boxplot$value, 
                                     color = "#9b87c4",
                                     fill = '#9b87c4')) +
  geom_violin(scale='area') +
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey", size=0.25)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_line(color='grey', size=0.25),
        legend.position = "none",
        text = element_text(size=7)) +
  scale_fill_manual(values=c("#9b87c4")) +
  scale_color_manual(values=c("#9b87c4"))

interface_a1bg_boxplot 


# 8.4 Supplementary Figure 5b,d - preparing data for abundance boxplots:

abundances_babad_boxplot <- data_babad[rowOfConditions_babad,1:3]
abundances_a1bg_boxplot <- data_a1bg[rowOfConditions_a1bg,1:3]
abundances_babag_boxplot <- data_babad[rowOfConditions_babag,1:3]


abundances_babad_boxplot <- pivot_longer(
  data = abundances_babad_boxplot,
  cols = colnames(abundances_babad_boxplot)[1:3],
  names_to = "interface",
)
abundances_a1bg_boxplot <- pivot_longer(
  data = abundances_a1bg_boxplot,
  cols = colnames(abundances_a1bg_boxplot)[1:3],
  names_to = "interface",
)

abundances_babag_boxplot <- pivot_longer(
  data = abundances_babag_boxplot,
  cols = colnames(abundances_babag_boxplot)[1:3],
  names_to = "interface",
)


# 8.5 Supplementary Figure 5b - abundance boxplots (a4b3d):

abundanceBoxplot_babad <- ggplot(data = abundances_babad_boxplot, 
                                 aes(x = abundances_babad_boxplot$interface, 
                                     y = abundances_babad_boxplot$value,
                                     fill = '#8cc67f',
                                     color = '#8cc67f')) +
  geom_boxplot(notch=TRUE, outlier.shape = NA) +
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey", size=0.25)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_line(color='grey', size=0.25),
        legend.position = "none",
        text = element_text(size=7)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#8cc67f")) +
  scale_color_manual(values=c("#8cc67f"))

abundanceBoxplot_babad



# 8.5 Supplementary Figure 5d - abundance boxplots (a1b3g):

abundanceBoxplot_a1bg <- ggplot(data = abundances_a1bg_boxplot, 
                                aes(x = abundances_a1bg_boxplot$interface, 
                                    y = abundances_a1bg_boxplot$value,
                                    fill = '#8cc67f',
                                    color = '#8cc67f')) +
  geom_boxplot(notch=TRUE, outlier.shape = NA) +
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey", size=0.25)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_line(color='grey', size=0.25),
        legend.position = "none",
        text = element_text(size=7)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#8cc67f")) +
  scale_color_manual(values=c("#8cc67f"))

abundanceBoxplot_a1bg



# 8.6 Mean abundances identified

# mean(subset(abundances_babad_boxplot, interface == 'alpha')$value)
# [1] 0.7646057
# mean(subset(abundances_babad_boxplot, interface == 'beta')$value)
# [1] 0.1224211
# mean(subset(abundances_babad_boxplot, interface == 'delta')$value)
# [1] 0.1129732
# 
# mean(subset(abundances_a1bg_boxplot, interface == 'alpha')$value)
# [1] 0.6183674
# mean(subset(abundances_a1bg_boxplot, interface == 'beta')$value)
# [1] 0.2105619
# mean(subset(abundances_a1bg_boxplot, interface == 'delta')$value)
# [1] 0.1710707
# 
# mean(subset(abundances_babag_boxplot, interface == 'alpha')$value)
# [1] 0.7901235
# mean(subset(abundances_babag_boxplot, interface == 'beta')$value)
# [1] 0.01234568
# mean(subset(abundances_babag_boxplot, interface == 'delta')$value)
# [1] 0.1975309


# 8.7 Saving boxplots and violin plots for Supplementary Fig. 5b-d

# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babad_interface_boxplot.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(interface_babad_boxplot)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/a1bg_interface_boxplot.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(interface_a1bg_boxplot)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/babad_abundance_boxplot.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(abundanceBoxplot_babad)
# dev.off()
# 
# pdf('/cephfs/asente/extrasynaptics/simulation-analysis/a1bg_abundance_boxplot.pdf', 
#     width = 1, height = 2, onefile = F) 
# print(abundanceBoxplot_a1bg)
# dev.off()


