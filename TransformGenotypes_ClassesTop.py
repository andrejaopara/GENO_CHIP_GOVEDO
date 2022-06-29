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
import commands
import tempfile
import pandas as pd
import sys


def remove_from_zip(zipfname, *filenames):
    """
    This is a function to remove files from a zipfile
    :param zipfname: name of the file to unzip
    :param filenames: files to remove
    :return: creates a new zip with the files removed
    """
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
date = raw_input("Vnesi datum [danes, brez presledkov]: ")
pasma = raw_input("Vnesi pasmo [Govedo(vse)/Rjava/Crnobela/Lisasta]: ")
AlleleFormat = raw_input("Vnesi način kodiranja alelov [top / forward / ab]: ")
Origin = raw_input("Vnesi vir [1-100]")
zip_file = raw_input("Vnesi ime zip datoteke z genotipi: ")
merge_ask = raw_input("Ali hočeš genotipe združiti s prejšnjimi genotipi (po čipu)? [Y/N] ")
sort_finalReport = raw_input("Ali hočeš sortirat FinalReport datoteko? [Y/N]")


# Ask what action does the user want to perform
parentageTest = raw_input("Ali hočeš izvleči SNP-e za preverjanje starševstva (razširjen set)?  [Y/N] ")
# Ask whether you want to remove original zip
rmOriginalZip = raw_input('Odstranim originalno zip datoteko? [Y/N] ')
# Create directory path to hold current temp genotype files within Genotipi_DATA and breed directory
tempDir = "/home/andreja/OBDELAVA_GENOTIPOV/Genotipi_DATA/" + pasma + "_TEMP/Genotipi_" + str(date) + "/"
# PEDDAROW directory
peddarow = "/home/andreja/OBDELAVA_GENOTIPOV/GENO_CHIP_GOVEDO"
# Zip latest
Zip_lat = "/home/andreja/OBDELAVA_GENOTIPOV/Genotipi_DATA/Genotipi_latest/" + pasma + "/Top/ZipGenoFiles/"
# Directory of the PLINK package directories and files
PLINKDIR = '/home/andreja/OBDELAVA_GENOTIPOV/Genotipi_DATA/Genotipi_latest/' + pasma + '/Top/'
# Path to the plink software
plinkSoftware = '~/bin/plink_linux_1.9_x86_64/plink' #naj bo pot/karkoliže/plink
# Path to Zanardi
ZanDir = "/home/andreja/OBDELAVA_GENOTIPOV/GENO_CHIP_GOVEDO/Zanardi/"
# Path to directory with essential files (usually the git directory)
CodeDir = "/home/andreja/OBDELAVA_GENOTIPOV/GENO_CHIP_GOVEDO"
# The directory of the downloaded genotype packages
DownloadDir = "/home/andreja/Downloads"


# File with IDs and seq for the animals
Breed_IDSeq = "/home/andreja/OBDELAVA_GENOTIPOV/GENO_CHIP_GOVEDO/" + pasma + "_seq_ID.csv"


# Name of the genotype zip file
zipPackage = zip_file
#########################################################################################################
##########################################################################################################
##########################################################################################################
# create dictionaries
##########################################################################################################
##########################################################################################################

# Create a dictionary of the number of SNPs and corresponding chip names
chips = GenFiles.chips

# Create empty dictionaries to hold the information
GenoFile = defaultdict(set)
SampleIDs = defaultdict(list)
AllInfo = []

# Dictionary to hold download date of the genotype package
DateDownloaded = defaultdict(list)
DateGenotyped = defaultdict(list)



# Read in animal ID / Seq / DateOfBirth / SexCode table
# Create a dictionary
Breed_IDSeq_Dict = defaultdict()
with open(Breed_IDSeq, 'rt') as IDSeq:
    reader = csv.reader(IDSeq, delimiter=',')
    for line in reader:
        Breed_IDSeq_Dict[line[0]] = line[1:]

############################################################################################################
#############################################################################################################
# Create a directory with the current date for temp genotype manipulation
if not os.path.exists(tempDir):
    os.makedirs(tempDir)

# Change current working directory to the created directory
os.chdir(tempDir)
# Copy the zip package into the "date" directory
shutil.copy(DownloadDir + "/" + zipPackage, tempDir)

# Create a onePackage object from the zip file
onePackage = GenFiles.genZipPackage(zipPackage)
# Extract the final report, SNPMap and SampleMap
onePackage.extractFinalReport()
onePackage.extractSNPMap()
onePackage.extractSampleMap()

# Sort the final report by animal if required
if sort_finalReport == "Y":
    os.system("bash " + CodeDir + "/SortFinalReport.sh " + onePackage.finalreportname + " " + CodeDir)

# Print the package statistics
print("Name of the package is: " + onePackage.name)
print("Name of the SNPmap is: " + onePackage.snpmapname)
print("Name of the Sample Map is: " + onePackage.samplemapname)
print("Name of the FinalReport is: " + onePackage.finalreportname)


# check for error IDs and replace the prior identified errouneous IDs
print("Obtaining spurious IDs")
# Read in the ErrorIDs file, in which we specify the erroneous IDs and their correct versions
replaceIDs = open(CodeDir + "/ErrorIDs_genotipi.txt").read().strip().split("\n")
replaceIDs = [tuple(replaceIDs[x].split(",")) for x in range(len(replaceIDs))]
# The extractErrorNames() function extracts the wrong IDs by pre-defined criteria --> check the function
#errorIDs = onePackage.extractErrorNames()  # extract Sample Names if they exist - they shouldnt be in the file
# to samo, če ti samo prav popravi!!!!!!!!!!!!!

# Create a list of wrong IDs
rID = [x for (x,y) in replaceIDs]
#errorIDs = [(x, y) for (x, y) in errorIDs if str(x) not in rID]
#errorIDs = errorIDs + replaceIDs
errorIDs = replaceIDs
print("Spurious IDs are: ")
print(errorIDs)

# Replace the wrong IDs in the final report file
print("Replacing spurious IDs")
if errorIDs:
    s = open(onePackage.samplemapname).read()
    sF = open(onePackage.finalreportname).read()

    for i in errorIDs:
        s = s.replace('\t' + str(i[0]), '\t' + str(i[1]))
        sF = sF.replace('\t' + str(i[0]), '\t' + str(i[1]))

    f = open(onePackage.samplemapname, 'w')
    f.write(s)
    f.close()
    fF = open(onePackage.finalreportname, 'w')
    fF.write(sF)
    fF.close()
    print("Successfully updated FinalReport and SampleMap.")

# Copy pedda.param and python script to the current directory
print("Preparing peddar.param parameter file")
shutil.copy((peddarow + "/peddar.param"), "peddar.param")
shutil.copy((peddarow + "/pedda_row.py"), "pedda_row.py")

# Replace strings in the param file for peddarow with shell command
os.system(
    'sed -i "s|test_FinalReport.txt|' + onePackage.finalreportname + '|g" peddar.param')  # insert FinalReport name into peddar.param
os.system(
    'sed -i "s|Dominant |Dominant_|g" ' + onePackage.finalreportname)  # problem Dominant Red with a space
os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.name + '_SNP_Map.txt')  ##problem Dominant Red with a space
os.system('sed -i "s/test_outputfile/"' + onePackage.name + '"/g" peddar.param')  # insert OutPut name into peddar.param
os.system(
    'sed -i "s/test_SNPMap.txt/"' + onePackage.name + '_SNP_Map.txt' + '"/g" peddar.param')  # insert SNPMap name into peddar.param
os.system(
    'sed -i "s/AlleleFormat/"' + AlleleFormat + '"/g" peddar.param')  # insert desired AlleleFormat name into peddar.param
os.system('sed -i "s/TEST/' + pasma + '/g" peddar.param')
os.system("python2.7 pedda_row.py")  # transform into ped and map file



# #ABFORMAT
# shutil.copy((peddarow+"/peddar.param"), "peddar.param")
# shutil.copy((peddarow+"/pedda_row.py"), "pedda_row.py")
# #replace strings with shell command
# os.system('sed -i "s|test_FinalReport.txt|'+ onePackage.name+"_FinalReport.txt" + '|g" peddar.param') #insert FinalReport name into peddar.param
# os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.name+"_FinalReport.txt") #problem Dominant Red with a space
# os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.name+'_SNP_Map.txt') ##problem Dominant Red with a space
# os.system('sed -i "s/test_outputfile/"'+onePackage.name+"_AB"+'"/g" peddar.param') #insert OutPut name into peddar.param
# os.system('sed -i "s/test_SNPMap.txt/"'+onePackage.name+'_SNP_Map.txt'+'"/g" peddar.param') #insert SNPMap name into peddar.param
# os.system('sed -i "s/AlleleFormat/"'+"ab"+'"/g" peddar.param') #insert desired AlleleFormat name into peddar.param
# os.system('sed -i "s/TEST/"'+pasma+'"/g" peddar.param')
# os.system("python2.7 pedda_row.py") #transform into ped and map file

# create a new zip file with corrected error names
# shutil.move(onePackage.name+'_Sample_Map.txt', 'Sample_Map.txt') #rename extracted SampleMap

# Create a new zip with corrected files
print("Adding correct files to .zip")
with zipfile.ZipFile(onePackage.name + '_FinalReport.zip', 'w', zipfile.ZIP_DEFLATED, allowZip64 = True) as myzip:
    myzip.write(onePackage.finalreportname)  # create new FinalReport zip
with zipfile.ZipFile(onePackage.name + '_Sample_Map.zip', 'w', zipfile.ZIP_DEFLATED, allowZip64 = True) as myzip:
    myzip.write(onePackage.samplemapname)  # create new Sample_Map.zip
remove_from_zip(onePackage.zipname, onePackage.oldfinalreportname)
remove_from_zip(onePackage.zipname, onePackage.oldsamplemapname)
with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED, allowZip64 = True) as z:
    z.write(onePackage.finalreportname)
with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED, allowZip64 = True) as z:
    z.write(onePackage.samplemapname)


# Create a pedFile object from the .ped file
pedfile = GenFiles.pedFile(onePackage.name + '.ped')
# Rename the .ped and .map
plinkfilename = onePackage.name.split("/")[-1] + "-" + str(len(pedfile.snps))
os.rename(onePackage.name + ".ped", plinkfilename + ".ped")
os.rename(onePackage.name + ".map", plinkfilename + ".map")
# Create a pedFile object from the .ped file and a mapFile object from the map file
pedfile = GenFiles.pedFile(plinkfilename + ".ped")
mapfile = GenFiles.mapFile(plinkfilename + ".map")


# Perform QC!
print("Peforming QC")
os.system("bash " + CodeDir + "/1_QC_FileArgs.sh " + pedfile.name + " " + pedfile.chip + " " + plinkSoftware + " " + tempDir)

if not os.path.isfile(pedfile.name + "_" + pedfile.chip + "_CleanIndsMarkers.ped"):
    print("No animals left after QC!")

# Extract the SNPs for parentage testing if required
# SNP_ISAG_196 =
# DOPOLNI, PROSIM
if parentageTest == 'Y':
    print("Extracting SNPs for parentage testing")
    pedFileQC = GenFiles.pedFile(pedfile.name + "_" + pedfile.chip + "_CleanInds.ped")
    pedFileQC.extractNamedSnpList("SNP_ISAG_196.txt", "196." + pedfile.chip + "-" + str(len(pedfile.snps)), plinkSoftware)
    pedFileQC.extractNamedSnpList("SNP_ICAR_554.txt", "554." + pedfile.chip + "-" + str(len(pedfile.snps)), plinkSoftware)
    pedFileQC.extractNamedSnpList("TRAIT_SNP_VERSAK", "TRAIT_SNP." + pedfile.chip + "-" + str(len(pedfile.snps)), plinkSoftware)
    pedFileQC.extractNamedSnpList("TRAIT_SNP_ICAR_554", "TRAIT_ICAR_554." + pedfile.chip + "-" + str(len(pedfile.snps)), plinkSoftware)

# Add file to the dictionary of chip files
GenoFile[pedfile.chip].add(pedfile.name)
DateDownloaded[date] += (pedfile.name)
DateGenotyped[onePackage.genodate] += [(x, pedfile.chip) for x in (pedfile.samples)]
AllInfo += [(x, pedfile.chip, pedfile.name, onePackage.genodate) for x in (pedfile.samples)]

notFound = []
for i in pedfile.samples:
    if i in Breed_IDSeq_Dict:
        SampleIDs[i] = [i, Breed_IDSeq_Dict.get(i)[0], onePackage.genodate, (pedfile.chip + "-" + str(len(pedfile.snps))), date, Origin]

    else:
        print("Sample ID " + i + " in " + pedfile.name + " not found!!!")
        notFound.append(i)


################################################################################################
###############################################################################################
# END OF THE LOOP
# Merge ped files if merge_ask = Y
# Create table for govedo
#############################################################################################
###############################################################################################
print("The number of genotyped animals is {}.".format(len(pedfile.samples)))
print("The number of found IDs in {}".format(len(SampleIDs)))
print("The number of genotype packages (different date of genotyping) is {}.".format(len(DateGenotyped)))
pd.DataFrame({"ID": list(set(pedfile.samples) ^ set(SampleIDs))}).to_csv("NotFoundIDs.csv", index=None)


# Create a table of individuals for govedo
# Columns are seq, chip, date genotyped
GenotypedInd = pd.DataFrame.from_dict(SampleIDs, orient='index', dtype=None)
GenotypedInd.columns = ['ID', 'ZIV_ID_SEQ', 'GenoDate', 'Chip', 'DownloadDate', 'Namen']
imiss = pd.read_csv(tempDir + pedfile.name + "_" + pedfile.chip + ".imiss", sep="\s+")[["IID", "F_MISS"]]
imiss.columns = ['ID', "F_MISS"]

Tabela = pd.merge(GenotypedInd, imiss, on="ID")
Tabela.to_csv(path_or_buf=tempDir + str(onePackage.genodate) + 'GovedoInd.csv', sep=",", index=False)

# Create a reduced table to import into the database
GenotypedInd_reduced = GenotypedInd[['ID', 'GenoDate', 'Chip', 'Namen']]
Tabela_reduced = pd.merge(GenotypedInd_reduced, imiss, on="ID")
Tabela_reduced = Tabela_reduced[['ID', 'GenoDate', 'Chip', 'F_MISS', 'Namen']]
Tabela_reduced.to_csv(path_or_buf=tempDir + str(onePackage.genodate) + 'GovedoInd_reduced.csv', sep=",", index=False)

print("Created table for Govedo: " + pasma)

os.system("rm *_FinalReport.txt *_FinalReport.zip")

# Merge the genotypes with the rest of the genotypes from the same chip
if merge_ask == "Y":
    # merge is outside the loop
    # merge all the chips needed updating
    chip = pedfile.chip
    pedFILE = tempDir + pedfile.pedname
    mapFILE = tempDir + mapfile.mapname
    if not os.path.exists(PLINKDIR + chip):
        os.makedirs(PLINKDIR + chip)
    shutil.copy(pedFILE, PLINKDIR + chip)
    shutil.copy(mapFILE, PLINKDIR + chip)
    os.chdir(PLINKDIR + chip)
    shutil.copy("/home/andreja/OBDELAVA_GENOTIPOV/GENO_CHIP_GOVEDO/PARAMFILE.txt", PLINKDIR + chip)
    if not os.path.isfile(PLINKDIR + chip + '/PLINK_MERGED.ped'):
        shutil.move(pedFILE, "PLINK_MERGED.ped")
        shutil.move(mapFILE, "PLINK_MERGED.map")
    if os.path.isfile(PLINKDIR + chip + '/PLINK_MERGED.ped'):
        mergeChipCommand = plinkSoftware + " --file PLINK_MERGED --cow --merge-list {0} --recode --out PLINK_MERGED".format(
            'MergeChip.txt')
        with open('MergeChip.txt', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=" ")
            writer.writerow([pedFILE, mapFILE])

    status, output = commands.getstatusoutput(mergeChipCommand)  # merge with plink

    if status == 0:
        print
        "Successfully merged " + chip + " " + PLINKDIR + " " + chip
    else:
        print
        "Merging went wrong, error: " + str(status)



