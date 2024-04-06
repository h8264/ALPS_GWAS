library(TwoSampleMR)
library(MRInstruments)
library(data.table)

exposure_data<- as.data.frame(fread("ALPS.txt"))
exposure_data<-format_data(exposure_data,
                           snp_col = "SNP",beta_col = "BETA",
                           se_col = "SE",effect_allele_col = "A1",
                           other_allele_col = "A2",
                           pval_col = "P",
                           eaf_col = "MAF",
                           chr_col = "CHR")
exposure_data<-exposure_data[exposure_data$pval.exposure<5e-8,]
exposure_data_clumped<-clump_data(exposure_data)
outcome_file<-paste("outcome_file.txt",file_name,sep = "")
outcome_data<-read_outcome_data(snps=exposure_data_clumped$SNP,
                                filename=outcome_file,
                                sep=" ",snp_col="SNP",
                                beta_col="BETA",se_col="SE",
                                effect_allele_col="A1",
                                other_allele_col="A2",
                                chr_col="CHR",
                                pval_col = "P")
dat<-harmonise_data(exposure_dat=exposure_data_clumped,outcome_dat = outcome_data)
results<-mr(dat)
dat.mr_heterogeneity<-mr_heterogeneity(dat)
dat.mr_pleiotropy_test<-mr_pleiotropy_test(dat)
