from pathlib import Path
from kmers import runKmers
from GeneCompare import geneCompare
from fileManager import *
import time
from progressBar import ProgressBar
from Stats import * 




testGene = getOrganismDir("testgene.fna")
gene1 = getOrganismDir("Ecoli.fna")
gene2 = getOrganismDir("VespaMandarinia.fna")
gene3 = getOrganismDir("AnophelesOryzalimnetes.fna")
gene4 = getOrganismDir("HapalochlaenaMaculosa.fna")



def main():
    initLog()
    
    
    gene = testGene
    
    runKmers(gene, 4, True, True) #select gene, k and the boolean parameters determining whitch k-mers algoriths are going to be used 1: freq, 2: generateall
    #geneCompare(gene, gene)
    #genome_len(gene)
    
 
    print("\n>> Programm finished.\n")


main()






