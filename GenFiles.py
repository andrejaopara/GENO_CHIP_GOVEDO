# -*- coding: utf-8 -*-
import os
from collections import defaultdict
import csv
import pandas as pd
import zipfile
import shutil
import re


chips = {19720: "GGPv02", 
26145: "GGPv03", 
26151: "GGPv03", 
30105: "GGPv04",
30106: "GGPv04", 
76883: "HD" , 
138892: "HDv02", 
139376: "HDv02", 
54001:"50Kv01" , 
54609: "50Kv02",
51274: "IDBv03",
52445: "IDBv03",
49629: "Versa50K",
49702: "Versa50K",
49706: "Versa50K",
49718: "Versa50K",
49740: "Versa50K",
6909: "IlluminaLD"
         }

TraitSNPs = {
'Caseins':['BCNAB', 'BCNAB_2', 'BCNAB_3','GNSC319','GNSC319_3','GNSC319_B1','GNSC355','GNSC355_3','GNSC355_B1','KappaCasein12951_1','KappaCasein12951_2','KappaCasein12951_3'],
'BetaCasein' : ['BCNAB', 'BCNAB_2', 'BCNAB_3'], 
'KappaCasein': ['GNSC319','GNSC319_3','GNSC319_B1','GNSC355','GNSC355_3','GNSC355_B1','KappaCasein12951_1','KappaCasein12951_2','KappaCasein12951_3'],
'Citrulinemia':[]}



"""
To je samo za ustvarit Sifrant SNPov,potem, ko ga imaš, ga skopiraj v datoteko
SNPSifrant="/home/jana/Genotipi/ParentalVerification_SNPSNP/Sifrant_SNP.csv"

SNPSifrant_Dict = defaultdict(set)
with open(SNPSifrant, 'rb') as SNP_Sifrant:
    SNP_Sifrant.readline() #skip the first line
    reader = csv.reader(SNP_Sifrant, delimiter=',')
    for line in reader:
        SNPSifrant_Dict[line[0]].add(line[1])
        
"""
SNP800Sifrant_Dict = {}

class genZipPackage:
    def __init__(self, zipDatoteka):
        self.zipFile=zipfile.ZipFile(zipDatoteka)
        self.zipname=zipDatoteka
        self.name=zipDatoteka.strip(".zip").strip(".ZIP")
        self.sernum=zipDatoteka.strip(".zip").strip('Matija_Rigler_')
        self.genodate=str([i for i in re.findall('\d+', self.zipname) if '2013' in i or '2014' in i or '2015' in i or '2017'
                       in i or '2018' in i or '2016' in i or '2019' in i][0])
        self.infiles=self.zipFile.namelist()
        self.finalreportname=[s for s in self.infiles if "final" in s.lower()][0] if len([s for s in self.infiles if "final" in s.lower()]) == 1 else [s for s in self.infiles if "final" in s.lower()]
        self.samplemapname = [s for s in self.infiles if "sample_map" in s.lower()][0] if len([s for s in self.infiles if "sample_map" in s.lower()]) == 1 else [s for s in self.infiles if "sample_map" in s.lower()]
        self.snpmapname = [s for s in self.infiles if "snp_map" in s.lower()][0] if len([s for s in self.infiles if "snp_map" in s.lower()]) == 1 else [s for s in self.infiles if "snp_map" in s.lower()]

    def unzip(self):
        self.zipFile.extractall()
    
    def extractFinalReport(self): #extracts all FinalReports and extracts them / if there is no FinalReport, it returns notice
        if self.finalreportname:
            if self.finalreportname.endswith('zip'):
                self.zipFile.extract(self.finalreportname)
                zipfile.ZipFile(self.finalreportname).extractall()
                os.remove(os.getcwd()+'/'+self.finalreportname)
                if "/" in self.finalreportname:
                    frn = self.finalreportname.split("/")[1]
                    shutil.move(frn.strip(".zip") + ".txt", self.name + '_FinalReport.txt')
                else:
                    shutil.move(self.finalreportname.strip(".zip") + ".txt", self.name + '_FinalReport.txt')

                self.oldfinalreportname = self.finalreportname
                self.finalreportname = self.name + '_FinalReport.txt'
            else:
                self.zipFile.extract(self.finalreportname)
                shutil.move(self.finalreportname, self.name + '_FinalReport.txt')
                self.oldfinalreportname = self.finalreportname
                self.finalreportname = self.name + '_FinalReport.txt'
        else:
            return 'No FinalReport infile in ' + self.name
        
    def extractSampleMap(self):
        if self.samplemapname:
            if self.samplemapname.endswith('zip'):
                self.zipFile.extract(self.samplemapname)
                zipfile.ZipFile(self.samplemapname).extractall()
                #os.remove(self.samplemapname)
                if "/" in self.samplemapname:
                    frn = self.samplemapname.split("/")[1]
                    shutil.move(frn.strip(".zip") + ".txt", self.name + '_Sample_Map.txt')
                else:
                    shutil.move(self.samplemapname.strip("zip") + "txt", self.name+'_Sample_Map.txt')
                self.oldsamplemapname = self.samplemapname
                self.samplemapname = self.name + '_Sample_Map.txt'
            else:
                self.zipFile.extract(self.samplemapname)
                shutil.move(self.samplemapname, self.name + '_Sample_Map.txt')
                self.oldsamplemapname = self.samplemapname
                self.samplemapname = self.name + '_Sample_Map.txt'
        else:
            return 'No Sample_Map.zip in {} infiles'.format(self.name)
        
    def extractSampleNames(self):
        if self.samplemapname:
            self.extractSampleMap()
            #os.system('sed -i "s|SI  |SI|g" ' + self.name+'_Sample_Map.txt') #remove double spacing
            #os.system('sed -i "s|SI |SI|g" ' + self.name+'_Sample_Map.txt') #remove space in sample IDs
            sampleTable = pd.read_table(self.name+'_Sample_Map.txt')
            return list(sampleTable['Name'])
        else:
            return 'No Sample_Map.zip in {} infiles'.format(self.name)

    
    def extractErrorNames(self): #this creates a list of tuples - errorID, correctedID
        if self.samplemapname:
            if not os.path.isfile(self.samplemapname):
                self.extractSampleMap()
            #os.system('sed -i "s|SI  |SI|g" ' + self.name+'_Sample_Map.txt') #remove double spacing
            #os.system('sed -i "s|SI |SI|g" ' + self.name+'_Sample_
            # Map.txt') #remove space in sample IDs
            sampleTable = pd.read_table(self.samplemapname, sep="\t")
            names=sampleTable['ID']
            errornames=[]
            for i in names:
                if str(i).isdigit() and len(str(i)) == 7 and str(i)[0] != 0:
                    errornames.append((i, "SI0" + str(i)))
                if str(i).isdigit() and len(str(i)) == 8:
                    errornames.append((int(i), 'SI'+str(i)))
                if str(i).find("  ") != -1:
                    a = str(i).replace("  ", " ")
                #
                #     if a.find(" ") != -1:
                #         errornames.append((i,( "".join(a.split(" ")[0:2]))))
                # if i.find("  ") == -1 and i.find(" ") != -1 and "SI" in i:
                #     errornames.append((i, "".join(i.split(" ")[:2])))
                    #
                if '(' in str(i) or ')' in str(i):
                    pass
            return errornames
        else:
            return ('No Sample_Map.zip in {} infiles'.format(self.name))
                    
                  
    def extractSNPMap(self):
        if self.snpmapname:
            if self.snpmapname.endswith('zip'):
                self.zipFile.extract(self.snpmapname)
                zipfile.ZipFile(self.snpmapname).extractall()
                #os.remove(self.snpmapname)
                if "/" in self.snpmapname:
                    frn = self.snpmapname.split("/")[1]
                    shutil.move(frn.strip(".zip") + ".txt", self.name + '_SNP_Map.txt')
                else:
                    shutil.move(self.snpmapname.strip("zip") + "txt", self.name + '_SNP_Map.txt')
            else:
                self.zipFile.extract(self.snpmapname)
                shutil.move(self.snpmapname, self.name + '_SNP_Map.txt')
        else:
            return 'No SNP_Map.zip in {} infiles'.format(self.name)

    def checkSubDir(self):
        if any(ime.endswith("/") for ime in self.infiles):
            return "Zip file contains multiple genotype packages, subdirectories will be created at file extraction"
            
    def zipSubDir(self, rmOriginalZip): #zips the 6 files and removes the directory
        if self.checkSubDir():
            self.subDirNames=filter(lambda x: x.endswith("/"), self.infiles)
            self.unzip()
            for i in self.subDirNames:
                shutil.make_archive(i, "zip",os.getcwd()+'/'+i)
                shutil.rmtree(os.getcwd()+'/'+i)
            
            if rmOriginalZip=='Y':
                os.remove(self.zipname)
                
                                     
                                                          
                                                                                                    
class pedFile:
    def __init__(self, pedDatoteka):
        self.pedname=pedDatoteka
        self.name=pedDatoteka.strip(".ped")
        self.pedContent=open(pedDatoteka).read().strip("\n").split("\n") #here don't read it in as panda table since it takes much longer
        self.samples=[line.split(" ")[1] for line in self.pedContent]
        self.mapContent=open(pedDatoteka.strip(".ped") + ".map").read().strip("\n").split("\n")
        try:
            self.snps=[line.split(" ")[1] for line in self.mapContent]
        except:
            self.snps=[line.split("\t")[1] for line in self.mapContent]
        try:
            self.chip=chips[(len(self.pedContent[0].split(" "))-6)/2]
        except: 
            self.chip = len(self.snps)
        #sort ped file by Individuals
        
    
    def extractSNP(self, SNP):
        os.system("plink --file " + self.name + " --cow --extract-snp "+ SNP + " --recode --out " + SNP)

    def extractNamedSnpList(self, SNPList):
        os.system("plink --file " + self.name + " --cow --extract " + SNPList + " --recode --out " +
                  self.pedname + "_" + SNPList.strip(".txt").strip(".csv"))

    def extractSNPList(self, SNPList):
        if "SNPList.txt" in os.listdir(os.getcwd()):
            overwrite=raw_input("Existing SNPList.txt file in the current working directory.\
 Do you want to overwrite and proceed? \n Ped and Map files will also be overwritten. [Y/N] ")
            if overwrite=='N':
                return None
            if overwrite == 'Y':
                pass
        with open("SNPList.txt", "w") as f:
            for snp in SNPList:
                f.write(snp + "\n")
        os.system("plink --file " + self.name + " --cow --extract SNPList.txt --recode --out SNPList")

    def extractTraitSNPs(self, Trait):
        SNPonChip = [x for x in TraitSNPs[Trait] if x in self.snps]
        if SNPonChip:
            print("\n"*3 +  "{0} SNPs found on the chip({1}): {2}".format(Trait, len(SNPonChip), SNPonChip))
            with open(Trait + ".txt", "w") as f:
                for snp in TraitSNPs[Trait]:
                    f.write(snp + "\n")
            os.system("plink --file " + self.name + " --cow --extract "+Trait+".txt --recode --out " + Trait)
        if not SNPonChip:
            return "No {} SNPs found on chip".format(Trait)
            
    def extractParentalSNPs(self, number):
        if number in (100,200,800):
            with open("ParentalSNP_" + str(number) + ".txt", "w") as f:
                for i in eval("SNP"+str(number)+"Sifrant_Dict"):
                    f.write(i + "\n")
            os.system("plink --file " + self.name + " --cow --extract ParentalSNP_800.txt --recode --out ParentalSNP800_"+self.sernum)
        else:
            return "Non-standard number of SNPs for parental verification"
  
    def individualSNPsDF(self, sampleID):
        position=[i for i,x in enumerate(self.samples) if x == sampleID][0]
        return  pd.DataFrame({self.pedContent[position].split(" ")[1] : self.pedContent[position].split(" ")[6:]})
        
    def individualSNPsList(self, sampleID):
        position=[i for i,x in enumerate(self.samples) if x == sampleID][0]
        return self.pedContent[position].split(" ")[6:]
        
        
        
    
    
    
class mapFile:
    def __init__(self, mapDatoteka):
        self.mapname=mapDatoteka
        self.name=mapDatoteka.strip(".map")
        self.sernum=mapDatoteka.strip(".map").strip("Matija_Rigler_")
        self.mapContent=pd.read_table(mapDatoteka, header=None, sep=" ")
        if len(self.mapContent.columns) != 4:
            self.mapContent=pd.read_table(mapDatoteka, header=None, sep="\t")
        self.mapContent.columns = ["chr", "snp", "start", "stop"]
        self.snps=self.mapContent["snp"]
        try:
            self.chip=chips[(len(self.mapContent))]       
        except:
            self.chip=len(self.snps)
        
    def chrSNPs(self, chromosome):
        printSNPs = raw_input("The number of SNPs on chromosome {0} is {1}. ".format(chromosome, len(self.mapContent[self.mapContent.chr==str(chromosome)]))+"Do you want to display all the SNPs? [Y/N] ")
        if printSNPs == 'Y':
            return self.mapContent[self.mapContent.chr==str(chromosome)]
        else:
            pass
            
            
    def posSNP(self, chromosome, pos):
        return self.mapContent[(self.mapContent.chr==str(chromosome)) & (self.mapContent.stop==pos)]
              
        
class mapFileMerged:
    def __init__(self, mapDatoteka):
        self.mapname=mapDatoteka
        self.name=mapDatoteka.strip(".map")
        self.sernum=mapDatoteka.strip(".map").strip("Matija_Rigler_")
        self.mapContent=pd.read_table(mapDatoteka, header=None, sep=" ", names=("chr", "snp", "start", "stop"))
        self.snps=self.mapContent["snp"]
        try:
            self.chip=chips[(len(self.mapContent))]       
        except:
            self.chip=len(self.snps)
        
    def chrSNPs(self, chromosome):
        printSNPs = raw_input("The number of SNPs on chromosome {0} is {1}. ".format(chromosome, len(self.mapContent[self.mapContent.chr==str(chromosome)]))+"Do you want to display all the SNPs? [Y/N] ")
        if printSNPs == 'Y':
            return self.mapContent[self.mapContent.chr==str(chromosome)]
        else:
            pass
            
            
    def posSNP(self, chromosome, pos):
        return self.mapContent[(self.mapContent.chr==str(chromosome)) & (self.mapContent.stop==pos)]
       
