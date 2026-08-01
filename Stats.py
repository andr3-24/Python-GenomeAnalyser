from pathlib import Path
from kmersV2 import runKmers
from log import *
from hash_utils import hashUse

#Get organisms folder
organism_folder = BASE_DIR / "Organisms"

def genome_len(gene):
    genepath = organism_folder / gene
    log("Calculating total genome length")
    
    """
    Calculates the total genome length in base pairs (bp).
    Parameters:
        filepath (Path): Path to FASTA file
    Returns:
        int: Total number of bases
    """

    total_length = 0

    with open(genepath, "r") as file:

        for line in file:

            line = line.strip()

            # Ignore FASTA headers
            if line.startswith(">"):
                continue

            total_length += len(line)
    log("Total genome length:", total_length)     
    print("> Total genome length:", total_length) 
    return total_length