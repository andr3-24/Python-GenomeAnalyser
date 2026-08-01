from pathlib import Path
from log import *




#Locate the gene files
def getOrganismDir(genename):
    
    #Get organisms folder
    organism_folder = BASE_DIR / "Organisms"
    
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