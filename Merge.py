# -*- coding: utf-8 -*-
# This is a script to add newly genotyped individuals and downloaded GeneSeek zip file (Final Reports)
# to the existing database of the latest genotypes

# Pipeline:
# 1)define chip, dictionary, dictionary to hold chip: genotype package: animal ids, dictionary to hold genotype package : download date
# 2)create temp directory within breed_TEMP/DownloadDate
# 3)create directory if not existin, unzip file
# 4)for each genotype package: unzip Final report and SNP_Map, change spurious strings within the files and adjust peddar.param file
# run peddar_row to transform FinalReports to PLINK and MAP formats
# write names to dictionaries


import os
import zipfile
import shutil
from collections import defaultdict
import csv
import GenFiles
import subprocess
import commands
import tempfile
import pandas as pd


def remove_from_zip(zipfname, *filenames):
    tempdir = tempfile.mkdtemp()
    try:
        tempname = os.path.join(tempdir, 'new.zip')
        with zipfile.ZipFile(zipfname, 'r') as zipread:
            with zipfile.ZipFile(tempname, 'w') as zipwrite:
                for item in zipread.infolist():
                    if item.filename not in filenames:
                        data = zipread.read(item.filename)
                        zipwrite.writestr(item, data)
        shutil.move(tempname, zipfname)
    finally:
        shutil.rmtree(tempdir)


########################################################
# set directories and file names
########################################################
# date='09082018'
# pasma="Rjava"
# AlleleFormat="top"
# zip_file="we_mr_19042018_IDB191.zip"
# merge_ask='N'
# Ask the user for the current date (date of download) and breed
date = raw_input("Enter the date (today): ")
pasma = raw_input("Enter the breed [Rjava/Crnobela/Lisasta]: ")
AlleleFormat = raw_input("Enter the desired allele coding [top / forward / ab]: ")
zip_file = raw_input("Enter the name of the downloaded zip file: ")
merge_ask = raw_input("Do you want to merge newly downloaded genotypes to the Latest Genotypes files (by chip)? [Y/N] ")

# ask what action does the user want to perform
# action = raw_input("Do you want to extract SNPs for parental verification  [Y/N] ")
action = 'N'
if action == 'Y':
    PVSNPs = input("How many SNPs would you like to use for parental verification? ")
# ask whether you want to remove original zip
rmOriginalZip = raw_input('Remove original zip? [Y/N] ')
# create directory path to hold current temp genotype files within Genotipi_DATA and breed directory
tempDir = "/home/jana/Genotipi/Genotipi_DATA/" + pasma + "_TEMP/Genotipi_" + str(date) + "/"
# PEDDAROW directory
peddarow = "/home/jana/Genotipi/TransformGeno/SNPchimpRepo/source_codes/PEDDA_ROW/"
# Zip latest
Zip_lat = "/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/" + pasma + "/Top/ZipGenoFiles/"
# Zip_lat="/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Zip/"
PLINKDIR = '/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/' + pasma + '/Top/'
# path to Zanardi
ZanDir = "/home/jana/Genotipi/TransformGeno/Zanardi/"
CodeDir = "/home/jana/Genotipi/TransformGeno/"
DownloadDir = "/home/jana/Downloads/"


# file with IDs and seq for the animals
Breed_IDSeq = "/home/jana/Genotipi/TransformGeno/" + pasma + "_seq_ID.csv"
Mesne_IDSeq = "/home/jana/Genotipi/TransformGeno/Mesne_seq_ID.csv"

os.chdir(tempDir)
# name of the file
zipPackage = zip_file
#########################################################################################################
##########################################################################################################
##########################################################################################################
# create dictionaries
##########################################################################################################
##########################################################################################################

# create a dictionary of the number of SNPs and corresponding chip names
chips = GenFiles.chips


GenoFile = defaultdict(set)
SampleIDs = defaultdict(list)
PedFiles = defaultdict(list)
MapFiles = defaultdict(list)
PedFilesQC = defaultdict(list)
MapFilesQC = defaultdict(list)
AllInfo = []
# dictionary to hold downlaod date of the genotype package
DateDownloaded = defaultdict(list)
DateGenotyped = defaultdict(list)





onePackage = GenFiles.genZipPackage(tempDir + "/" + zipPackage)




# make pedfile a GenFiles pedFile object
# make pedfile a GenFiles pedFile object
pedfile = GenFiles.pedFile(onePackage.name + '.ped')
pedfilename = onePackage.name.split("/")[-1]
mapfile = GenFiles.mapFile(onePackage.name + '.map')

# Perform QC!
print("Peforming QC")
#os.system("cp " + CodeDir + "/1_QC_FileArgs.sh " + tempDir)
#os.system("bash ./1_QC_FileArgs.sh " + pedfilename + " " + pedfile.chip)
#PedFilesQC[pedfile.chip].append(tempDir + pedfilename + "_" + pedfile.chip + "_CleanIndsMarkers.ped")
#MapFilesQC[pedfile.chip].append(tempDir + pedfilename + "_" + pedfile.chip + "_CleanIndsMarkers.map")

# add file to the dictionary of chip files
PedFiles[pedfile.chip].append(pedfile.pedname)
MapFiles[pedfile.chip].append(mapfile.mapname)
GenoFile[pedfile.chip].add(pedfile.name)
DateDownloaded[date] += (pedfile.name)
DateGenotyped[onePackage.genodate] += [(x, pedfile.chip) for x in (pedfile.samples)]
AllInfo += [(x, pedfile.chip, pedfile.name, onePackage.genodate) for x in (pedfile.samples)]



if merge_ask == "Y":
    # merge is outside the loop
    # merge all the chips needed updating

    for i in PedFiles:
        if not os.path.exists(PLINKDIR + str(i)):
            os.makedirs(PLINKDIR + str(i))
        for pedfile, mapfile in zip(PedFiles[i], MapFiles[i]):
            shutil.copy(pedfile, PLINKDIR + str(i))
            shutil.copy(mapfile, PLINKDIR + str(i))
        os.chdir(PLINKDIR + str(i))
        shutil.copy("/home/jana/Genotipi/Genotipi_CODES/PARAMFILE.txt", PLINKDIR + i)
        pedToMerge = ",".join(PedFiles[i]).strip("'")
        mapToMerge = ",".join(MapFiles[i]).strip("'")
        if not os.path.isfile(PLINKDIR + i + '/PLINK_MERGED.ped'):
            mergeChipCommand = "plink --file {0} --cow --merge-list {1} --recode --out PLINK_MERGED".format(
                (PedFiles[i][0].strip(".ped")), 'MergeChip.txt')
            with open('MergeChip.txt', 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=" ")
                [writer.writerow(r) for r in
                 zip(PedFiles[i][1:], MapFiles[i][1:])]  # leave the first one out - that goes in the plink command line
        if os.path.isfile(PLINKDIR + i + '/PLINK_MERGED.ped'):
            mergeChipCommand = "plink --file PLINK_MERGED --cow --merge-list {0} --recode --out PLINK_MERGED".format(
                'MergeChip.txt')
            with open('MergeChip.txt', 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=" ")
                [writer.writerow(r) for r in zip(PedFiles[i], MapFiles[i])]

        status, output = commands.getstatusoutput(mergeChipCommand)  # merge with plink

        if status == 0:
            print
            "Successfully merged " + str(i) + " " + PLINKDIR + " " + i
        else:
            print
            "Merging went wrong, error: " + str(status)

for chip in PedFiles:
    PedFiles[chip] = [i.replace("ZipGenoFiles", "ZipGenoFiles/") for i in PedFiles[chip]]

for chip in MapFiles:
    MapFiles[chip] = [i.replace("ZipGenoFiles", "ZipGenoFiles/") for i in MapFiles[chip]]

# MERGE FOR QC-ed data!!!!

# for i in PedFiles:
#    if not os.path.exists(PLINKDIR+str(i)):
#        os.makedirs(PLINKDIR+str(i))
#    for pedfile, mapfile in zip (PedFilesQC[i], MapFilesQC[i]):
#        shutil.copy(pedfile, PLINKDIR+str(i))
#        shutil.copy(mapfile, PLINKDIR+str(i))
#    os.chdir(PLINKDIR+str(i))
#    shutil.copy("/home/jana/Genotipi/Genotipi_CODES/PARAMFILE.txt", PLINKDIR+i)
#    pedToMerge = ",".join(PedFilesQC[i]).strip("'")
#    mapToMerge = ",".join(MapFilesQC[i]).strip("'")
#    if not os.path.isfile(PLINKDIR+i+'/PLINK_MERGED_' + i + '_CleanIndsMarkers.ped'):
#        mergeChipCommand = "plink --file {0} --cow --merge-list {1} --recode --out {2}".format((PedFilesQC[i][0].strip(".ped")), 'MergeChip.txt', "PLINK_MERGED_" + i + "_CleanIndsMarkers")
#        with open('MergeChip.txt', 'w') as csvfile:
#            writer = csv.writer(csvfile, delimiter=" ")
#            [writer.writerow(r) for r in zip(PedFilesQC[i][1:], MapFilesQC[i][1:])] #leave the first one out - that goes in the plink command line
#    if os.path.isfile(PLINKDIR+i+'/PLINK_MERGED_' + i + '_CleanIndsMarkers.ped'):
#        mergeChipCommand = 'plink --file PLINK_MERGED_{0}_CleanIndsMarkers --cow --merge-list {1} --recode --out PLINK_MERGED_{0}_CleanIndsMarkers'.format(i, 'MergeChip.txt')
#        with open('MergeChip.txt', 'w') as csvfile:
#            writer = csv.writer(csvfile, delimiter=" ")
#            [writer.writerow(r) for r in zip(PedFilesQC[i], MapFilesQC[i])]
#
#    status, output = commands.getstatusoutput(mergeChipCommand) #merge with plink
#
#    if status == 0:
#        print "Successfully merged " + str(i) + " " + PLINKDIR + " " + i + "_CleanIndsMarkers"
#    else:
#       print "Merging went wrong, error: " + str(status)