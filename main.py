from pathlib import Path
from kmers import runKmers
from log import *
from GeneCompare import geneCompare
from hash_utils import hashUse
from Stats import * 
from fileManager import *





def importFile():
    #User imports the .fna file in the program 
    #Programm automatically creates the corresponding directory in the 'Organisms' folder
    #The new directory name is filename - the file extension 
    #e.g example.fna : Organisms\example\example.fna
    pass



testGene = getOrganismDir("testgene.fna")
gene1 = getOrganismDir("Ecoli.fna")
gene2 = getOrganismDir("VespaMandarinia.fna")
gene3 = getOrganismDir("AnophelesOryzalimnetes.fna")
gene4 = getOrganismDir("HapalochlaenaMaculosa.fna")



def main():
    gene = testGene
    
    initLog()
    runKmers(gene, 3, True, True) #select gene, k and the boolean parameters determining whitch k-mers algoriths are going to be used 1: freq, 2: generateall
    #geneCompare(gene, gene)
    #genome_len(gene)
    
    
    
    print("\n>> Programm finished.\n")


main()






