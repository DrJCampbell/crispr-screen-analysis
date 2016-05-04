#
#
#

setwd("/Users/jamesc/Dropbox/212_Jess_CRISPR_screens/")


# Sample 18 is actually sample 17 - check with Kerry/Ioannis
sample18 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0001.IonXpress_018.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)



sample19 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0002.IonXpress_019.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample20 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0003.IonXpress_020.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample21 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0004.IonXpress_021.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample22 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0005.IonXpress_022.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample23 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0006.IonXpress_023.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample24 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0007.IonXpress_024.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample25 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0008.IonXpress_025.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample26 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0009.IonXpress_026.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample27 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0010.IonXpress_027.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample28 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0011.IonXpress_028.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample29 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0012.IonXpress_029.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample29 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0012.IonXpress_029.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

sample30 <- read.table(
	file="R_2016_04_12_16_01_38_user_ICRProton-75-C81-2.C81-0013.IonXpress_030.fastq.adaptrm1.results.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)



# Get rid of controls (Olfr and PLK1)
# as makes no sense to consider with 
# samples. We find that all controls
# with a common sequence collapse to
# a single id (regardless of the plate
# they were on)

sample18_no_controls <- remove_controls(sample18)
sample19_no_controls <- remove_controls(sample19)
sample20_no_controls <- remove_controls(sample20)
sample21_no_controls <- remove_controls(sample21)
sample22_no_controls <- remove_controls(sample22)
sample23_no_controls <- remove_controls(sample23)
sample24_no_controls <- remove_controls(sample24)
sample25_no_controls <- remove_controls(sample25)


# ========================== #
#
# Negative selection screens
#
# ========================== #



# samples 21:25 are negative selection screens.
# 21 is t0,
# 22,23 are T=1 and T=2
# 24,25 are T=1 and T=2 plus BMN

sample21_pptm <- get_pptm(
	sample21_no_controls
	)
sample22_pptm <- get_pptm(
	sample22_no_controls
	)
sample23_pptm <- get_pptm(
	sample23_no_controls
	)
sample24_pptm <- get_pptm(
	sample24_no_controls
	)
sample25_pptm <- get_pptm(
	sample25_no_controls
	)

#
#
#

# 
# Z-scores without excluding low counts
# viabilty effect at 1wk
znorm_22_21_all <- znorm(
	x0=sample21_pptm,
	x1=sample22_pptm,
	min_counts=0
	)

# viabilty effect at 2wks
znorm_23_21_all <- znorm(
	x0=sample21_pptm,
	x1=sample23_pptm,
	min_counts=0
	)


# drug effect at 1wk
znorm_24_22_all <- znorm(
	x0=sample22_pptm,
	x1=sample24_pptm,
	min_counts=0
	)


# drug effect at 2wks
znorm_25_23_all <- znorm(
	x0=sample23_pptm,
	x1=sample25_pptm,
	min_counts=0
	)


#
# Z-scores after removing low count
# rows the t0/no-drug sample
#

# viabilty effect at 1wk
znorm_22_21 <- znorm(
	x0=sample21_pptm,
	x1=sample22_pptm
	)

# viabilty effect at 2wks
znorm_23_21 <- znorm(
	x0=sample21_pptm,
	x1=sample23_pptm
	)


# drug effect at 1wk
znorm_24_22 <- znorm(
	x0=sample22_pptm,
	x1=sample24_pptm
	)


# drug effect at 2wks
znorm_25_23 <- znorm(
	x0=sample23_pptm,
	x1=sample25_pptm
	)


#
# Histograms of Z-scores (as example)
# for two week drug and viability effect
#
pdf("Histograms_of_DE_VE_Scores_with_all_and_gt50_count_data.pdf", width=5, height=5)
# viabilty effect at 2wks
hist(
	znorm_23_21_all$x1_x0_zscore,
	xlab="VE Z-score (2 weeks)",
	main="All guides"
	)

# drug effect at 2wks
hist(
	znorm_25_23_all$x1_x0_zscore,
	xlab="DE Z-score (2 weeks)",
	main="All guides"
	)

# viabilty effect at 2wks
hist(
	znorm_23_21$x1_x0_zscore,
	xlab="VE Z-score (2 weeks)",
	main="Guides > 50 counts"
	)

# drug effect at 2wks
hist(
	znorm_25_23$x1_x0_zscore,
	xlab="DE Z-score (2 weeks)",
	main="Guides > 50 counts"
	)
dev.off()


#
# join up Z-scores with original data
# and write out
#

# viability at 1 wk
write.table(
	znorm_22_21,
	file="Sample22_21_filtered_viability_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# viabilty effect at 2wks
write.table(
	znorm_23_21,
	file="Sample23_21_filtered_viability_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# drug effect at 1wk
write.table(
	znorm_24_22,
	file="Sample24_22_filtered_drug_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# drug effect at 2wks
write.table(
	znorm_25_23,
	file="Sample25_23_filtered_drug_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


#
# write out 'all' zscores tables
#

#
# join up Z-scores with original data
# and write out
#

# viability at 1 wk
write.table(
	znorm_22_21_all,
	file="Sample22_21_all_viability_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# viabilty effect at 2wks
write.table(
	znorm_23_21_all,
	file="Sample23_21_all_viability_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# drug effect at 1wk
write.table(
	znorm_24_22_all,
	file="Sample24_22_all_drug_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# drug effect at 2wks
write.table(
	znorm_25_23_all,
	file="Sample25_23_all_drug_effect_Zscores_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


# ========================== #
#
# Positive selection screens
#
# ========================== #


sample19_no_controls <- remove_controls(sample19)
sample20_no_controls <- remove_controls(sample20)

sample26_no_controls <- remove_controls(sample26)
sample27_no_controls <- remove_controls(sample27)
sample28_no_controls <- remove_controls(sample28)
sample29_no_controls <- remove_controls(sample29)
sample30_no_controls <- remove_controls(sample30)

sample19_pptm <- get_pptm(
	sample19_no_controls
	)
sample20_pptm <- get_pptm(
	sample20_no_controls
	)

sample26_pptm <- get_pptm(
	sample26_no_controls
	)
sample27_pptm <- get_pptm(
	sample27_no_controls
	)
sample28_pptm <- get_pptm(
	sample28_no_controls
	)
sample29_pptm <- get_pptm(
	sample29_no_controls
	)
sample30_pptm <- get_pptm(
	sample30_no_controls
	)


#
# 26 vs 19
# SUM149 BMN positive selection (with PARP1 guides)
#

make_scatter_plot(
	sample19_pptm,
	sample26_pptm,
	filename="library_vs_sum149BMN.pdf",
	gene="PARP1",
	x0name="library",
	x1name="SUM149 BMN + PARP guides",
	main="PPTM counts"
	)

make_scatter_plot(
	sample19_pptm,
	sample26_pptm,
	filename="library_vs_sum149BMN_total_hits.pdf",
	gene="PARP1",
	response_colname="total.hits",
	x0name="library",
	x1name="SUM149 BMN + PARP guides",
	main="total counts"
	)

lm1_totalhits <- lm(
	log10(sample26_pptm$total.hits+0.5) ~ 	log10(sample19_pptm$total.hits+0.5)
	)

sample_19_26_with_residuals_lm1_totalhits <- cbind(
	sample19_pptm,
	sample26_pptm, 
	lm1_totalhits$residuals
	)

write.table(
	sample_19_26_with_residuals_lm1_totalhits,
	file="sample_19_26_with_residuals_total_hits_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


#
# 27 vs 19
# TOV21G 100nm VX970 positive selection (with PARP1 guides)
#

make_scatter_plot(
	sample19_pptm,
	sample27_pptm,
	filename="library_vs_TOV21G_100nm_VX.pdf",
	gene="no_gene",
	x0name="library",
	x1name="TOV21G 100nM VX",
	main="PPTM counts"
	)

make_scatter_plot(
	sample19_pptm,
	sample27_pptm,
	filename="library_vs_TOV21G_100nm_VX_total_hits.pdf",
	gene="no_gene",
	response_colname="total.hits",
	x0name="library",
	x1name="TOV21G 100nM VX",
	main="total counts"
	)

lm_TOV21G_100nm_VX_totalhits <- lm(
	log10(sample27_pptm$total.hits+0.5) ~ 	log10(sample19_pptm$total.hits+0.5)
	)

sample_19_27_with_residuals_lm1_totalhits <- cbind(
	sample19_pptm,
	sample26_pptm, 
	lm_TOV21G_100nm_VX_totalhits$residuals
	)

write.table(
	sample_19_27_with_residuals_lm1_totalhits,
	file="sample_19_27_with_residuals_lm1_TOV21G_100nm_VX_totalhits_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)

#
# 29 vs 19
# TOV21G 500nm VX970 positive selection (with PARP1 guides)
#

make_scatter_plot(
	sample19_pptm,
	sample29_pptm,
	filename="library_vs_TOV21G_500nm_VX.pdf",
	gene="no_gene",
	x0name="library",
	x1name="TOV21G 500nM VX",
	main="PPTM counts"
	)

make_scatter_plot(
	sample19_pptm,
	sample29_pptm,
	filename="library_vs_TOV21G_500nm_VX_total_hits.pdf",
	gene="no_gene",
	response_colname="total.hits",
	x0name="library",
	x1name="TOV21G 500nM VX",
	main="total counts"
	)

lm_TOV21G_500nm_VX_totalhits <- lm(
	log10(sample29_pptm$total.hits+0.5) ~ 	log10(sample19_pptm$total.hits+0.5)
	)

sample_19_29_with_residuals_lm1_totalhits <- cbind(
	sample19_pptm,
	sample29_pptm, 
	lm_TOV21G_500nm_VX_totalhits$residuals
	)

write.table(
	sample_19_29_with_residuals_lm1_totalhits,
	file="sample_19_29_with_residuals_lm1_TOV21G_500nm_VX_totalhits_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)


#
# 30 vs 19
# TOV21G 200nm 0.5m cells VX970 positive selection (with PARP1 guides)
#

make_scatter_plot(
	sample19_pptm,
	sample30_pptm,
	filename="library_vs_TOV21G_200nm_VX.pdf",
	gene="no_gene",
	x0name="library",
	x1name="TOV21G 200nM VX",
	main="PPTM counts"
	)

make_scatter_plot(
	sample19_pptm,
	sample30_pptm,
	filename="library_vs_TOV21G_200nm_VX_total_hits.pdf",
	gene="no_gene",
	response_colname="total.hits",
	x0name="library",
	x1name="TOV21G 200nM VX",
	main="total counts"
	)

lm_TOV21G_200nm_VX_totalhits <- lm(
	log10(sample30_pptm$total.hits+0.5) ~ 	log10(sample19_pptm$total.hits+0.5)
	)

sample_19_30_with_residuals_lm1_totalhits <- cbind(
	sample19_pptm,
	sample30_pptm, 
	lm_TOV21G_200nm_VX_totalhits$residuals
	)

write.table(
	sample_19_30_with_residuals_lm1_totalhits,
	file="sample_19_30_with_residuals_lm1_TOV21G_200nm_VX_totalhits_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)



#
# Compare sample 18 (without PARP1 guides???)
# to sample 19 (with PARP1 guides)
#

#
# 20 is the library with no PARP1 guides
# 18 is SUM149 BMN screen without PAPP1 guides
#

sample18_no_controls <- remove_controls(sample18)

sample18_pptm <- get_pptm(
	sample18_no_controls
	)

make_scatter_plot(
	sample20_pptm,
	sample18_pptm,
	filename="library_vs_SUM149_noPARP.pdf",
	gene="PARP",
	x0name="library (no PARP guides)",
	x1name="SUM149 BMN (no PARP guides)",
	main="PPTM counts"	
	)

make_scatter_plot(
	sample20_pptm,
	sample18_pptm,
	filename="library_vs_SUM149_noPARP_total_hits.pdf",
	gene="PARP",
	response_colname="total.hits",
	x0name="library (no PARP guides)",
	x1name="SUM149 BMN (no PARP guides)",
	main="total counts"
	)

lm_SUM149_noPARP_totalhits <- lm(
	log10(sample18_pptm$total.hits+0.5) ~ 	log10(sample20_pptm$total.hits+0.5)
	)

sample_18_20_with_residuals_lm1_totalhits <- cbind(
	sample18_pptm,
	sample20_pptm, 
	lm_SUM149_noPARP_totalhits$residuals
	)

write.table(
	sample_18_20_with_residuals_lm1_totalhits,
	file="sample_18_20_with_residuals_lm1_SUM149_noPARP_totalhits_160504.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t",
	quote=FALSE
	)

