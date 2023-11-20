#########################
# Transform final report to PLINK files
#########################
# Note: This is *one* of the possible solutions
#       Any other way could be used that creates the same file structure
#########################

# Clear workspace
rm(list = ls())
# Set working directory
setwd("C:/Users/jernejb/Desktop/snp_data")
#load packages
library(tidyverse)
library(dplyr)

# read in final report file

#
# tidyverse - keeps spaces in column names and these are between single quotation marks 'column Name'  
finalReportslo <-  read_delim("Aslo_2022_1-2_FinalReport.txt.gz", delim = "\t", skip = 9, col_names = T)
# if ncol(finalReport) = 0 set delim parameter to , (comma)

###############################
# Your goal is to create lgen, fam and map files
# Some information might be missing in the final report, so you need to replace them 
###############################
# Scenario 1 - You have all the info in the Final Report - Jackpot!
###############################
# Fam file
finalReportslo %>%
  distinct(`Sample Name`) %>%
  mutate(FID = "AMEL", sire = 0, dam = 0, sex = 0, phenotype = -9) %>%
  relocate(`Sample Name`, .after = FID) %>%
  write_delim("Aslo.fam", col_names = F)

# Lgen file
finalReportslo %>%
  mutate(FID = "AMEL") %>%
  select(FID, `Sample Name`, `SNP Name`, `Allele1 - AB`, `Allele2 - AB`) %>%
  write_delim("Aslo.lgen", col_names = F)

# Map file
finalReportslo %>%
  distinct(`SNP Name`, .keep_all = TRUE) %>%
  mutate(morgan = 0) %>%
  dplyr::select(Chr, `SNP Name`, morgan, Position) %>%
  write_delim("Aslo.map", col_names = F)

# change to ped file with PLINK 
system("plink --chr-set 16 --nonfounders --allow-no-sex --lfile Aslo --missing-genotype - --output-missing-genotype 0 --recode --out Aslo")

##############################
#CRO SAMPLES
##############################
finalReportCro <-  read_delim("Acro_2022_1-3_FinalReport.txt.gz", delim = "\t", skip = 9, col_names = T)

# Fam file
finalReportCro %>%
  distinct(`Sample Name`) %>%
  mutate(FID = "AMEL", sire = 0, dam = 0, sex = 0, phenotype = -9) %>%
  relocate(`Sample Name`, .after = FID) %>%
  write_delim("Acro.fam", col_names = F)

# Lgen file
finalReportCro %>%
  mutate(FID = "AMEL") %>%
  select(FID, `Sample Name`, `SNP Name`, `Allele1 - AB`, `Allele2 - AB`) %>%
  write_delim("Acro.lgen", col_names = F)

# Map file
finalReportCro %>%
  distinct(`SNP Name`, .keep_all = TRUE) %>%
  mutate(morgan = 0) %>%
  dplyr::select(Chr, `SNP Name`, morgan, Position) %>%
  write_delim("Acro.map", col_names = F)

# change to ped file with PLINK 
system("plink --chr-set 16 --nonfounders --allow-no-sex --lfile Acro --missing-genotype - --output-missing-genotype 0 --recode --out Acro")


system(paste0("plink --file Aslo --chr-set 16  --merge Acro",
              " --recode --out BeePop"))


