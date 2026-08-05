from pathlib import Path
from log import *
from hash_utils import hashUse

#Get organisms folder
organism_folder = BASE_DIR / "Organisms"

def genomeLen(gene):
    log("Calculating total genome length")
    genepath = organism_folder / gene
    tLen = 0

    """
    Calculates the total genome length in base pairs (bp).
    Parameters:
        filepath (Path): Path to FASTA file
    Returns:
        int: Total number of bases
    """
    
    with open(genepath, "r") as file:

        for line in file:

            line = line.strip()

            # Ignore FASTA headers
            if line.startswith(">"):
                continue

            tLen += len(line)
    log("Total genome length:", tLen)     
    return tLen

def genomeTotalLines(gene):
    genepath = organism_folder / gene
    
    log("Calculating total genome lines")

    totalines = 0
    with open(genepath, "r") as file:

        for line in file:
            totalines+=1

    log("Total genome lines:", totalines)      
    return totalines


def inputFileSize(gene):
    """
    Returns the size of a file in bytes.
    """
    
    genepath = organism_folder / gene
    
    if genepath.exists():
        return genepath.stat().st_size 
    
    return None