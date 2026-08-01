from pathlib import Path
from kmers import runKmers
from log import *
from GeneCompare import geneCompare
from hash_utils import hashUse
from Stats import * 



#Get organisms folder
organism_folder = BASE_DIR / "Organisms"

def importFile():
    #User imports the .fna file in the program 
    #Programm automatically creates the corresponding directory in the 'Organisms' folder
    #The new directory name is filename - the file extension 
    #e.g example.fna : Organisms\example\example.fna
    pass

#Locate the gene files
def getFilePath(genename):
    
    #Each .fna file exist in its own directory. In the 'genefolder' we track the path for each .fna file 
    #.stem -> cuts the file extension for each file
    
    genefolder = organism_folder / Path(genename).stem 
    
    #finally, the file_path variable holds the filename path
    
    file_path = genefolder / genename

    if not file_path.exists():
        raise FileNotFoundError(
            f"{genename} was not found in {organism_folder}"
        )

    return file_path

testGene = getFilePath("testgene.fna")
gene1 = getFilePath("Ecoli.fna")
gene2 = getFilePath("VespaMandarinia.fna")
gene3 = getFilePath("AnophelesOryzalimnetes.fna")
gene4 = getFilePath("HapalochlaenaMaculosa.fna")



def main():
    initLog()
    runKmers(testGene, 3, False, False) #select gene, k and the boolean parameters determining whitch k-mers algoriths are going to be used
    #geneCompare(gene1, gene2)
    #genome_len(gene1)
    
    
    
    print("Programm Finished.")


main()






