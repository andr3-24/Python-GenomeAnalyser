import os
from collections import Counter
from pathlib import Path
from log import log
from hash_utils import hashUse
from fileManager import *


'''
The k-mer algorithm extracts all possible subsequences of length k from a DNA sequence.
Each k-mer represents a short fragment of the genome and can be used for sequence
comparison, similarity analysis, and genome characterization.
'''

def getKmers(inputseq,k,everykmer): 
    kmercounter=0
    if ((os.path.getsize(inputseq) != 0) and (os.path.getsize(inputseq) >= k)): #If file is not empty (((os.path.getsize(inputseq) != 0)), and if it is greater than k, proceed,
        prev = ""
        with open(inputseq, "r") as inputfile, open(everykmer, "w") as outfile:
            for line in inputfile:
                if line.startswith(">"):
                    continue      
                line = line.strip().upper()
                fulline = prev + line
                
                for i in range (0, len(fulline)-k+1,1):
                    outfile.write(fulline[i:i+k] + "\n")
                    kmercounter+=1
                    # Instead of using an embedded loop, use fulline[i:i+k]. 
                    # It reads the DNA bases from index i to i+k-1, which is a significantly faster method.
                        #kmer = fulline[i:i+k]
                        #kmers.append(kmer)
                prev = fulline[-(k-1):]  
        hashUse(everykmer, 0)     #create sha256 file to secure authenticity
        log("K-mers file is ready.")
        #print(">> K-mers file is ready.")
        log(kmercounter, "k-mers of size", k , "have been found.")
        #print(kmercounter, "k-mers of size", k , "have been found.")
    else:
        log("Error while trying to read the file.")
        #print("Error while trying to read the file.")

def countkmerfrequency(inputseq,k,kmerfr):
            """
            Reads a genome file and creates a k-mer frequency file.

            The function:
            1) Reads the genome directly without creating an intermediate k-mer file (unless this action is needed -> getKmers function).
            2) Validates nucleotide symbols.
            3) Counts k-mer occurrences.
            4) Creates an output file containing k-mers sorted by frequency.
            """

            valid_bases = set("ATCGN")
            counts = Counter() # Creates an empty Counter object, which is a specialized dictionary
                               # from the collections module used for counting occurrences of elements.
                               
            ambiguous_bases = 0  # Counts the total number of ambiguous nucleotide symbols (e.g. 'N') found in the input genome.                   
            skipped_kmers = 0    # Counts the number of generated k-mers that were skipped because
                                 # they contain one or more ambiguous nucleotide symbols.
                                 
            # Check if file exists and is large enough
            if not os.path.exists(inputseq):
                raise FileNotFoundError(f"Genome file not found: {inputseq}") # Raises an exception and stops execution if the input genome file does not exist.
                log("Genome file not found: {inputseq}")
            if (os.path.getsize(inputseq) < k) or (os.path.getsize(inputseq) == 0):
                raise ValueError("Genome file is empty or smaller than given k.") # Raises an exception and stops execution if the input genome file is empty.
                log("Genome file is empty or smaller than given k.")
            
            # Stores the last (k-1) bases from the previous sequence line
            prev = ""
            
            with open(inputseq, "r") as inputfile:

                for line in inputfile:

                    # Ignore FASTA headers
                    if line.startswith(">"):
                        continue

                    sequence = line.strip().upper()
                    
                    # Skip empty lines
                    if not sequence:
                        continue

                    # Check for invalid nucleotide symbols
                    invalid = set(sequence) - valid_bases

                    if invalid:
                        raise ValueError(
                            f"Invalid nucleotide symbols detected: {invalid}"
                        )
                        log(f"Invalid nucleotide symbols detected: {invalid}")

                    # Count ambiguous bases
                    ambiguous_bases += sequence.count("N")
                    
                    # Join the current line with the last (k-1) bases of the previous line
                    full_sequence = prev + sequence
                    
                    # Generate k-mers directly
                    for i in range(len(full_sequence) - k + 1):
            
                        kmer = full_sequence[i:i+k]
            
                        # Ignore k-mers containing unknown bases
                        if "N" in kmer:
                            skipped_kmers += 1
                            continue
            
                        counts[kmer] += 1
                        
                        # Store the last (k-1) bases for the next iteration
                    if len(full_sequence) >= k - 1:
                        prev = full_sequence[-(k-1):]
                    else:
                        prev = full_sequence

            if not counts:
                raise ValueError(
                    "No valid k-mers were generated. Check genome length or k value."
                )
                log("No valid k-mers were generated. Check genome length or k value.")

            # Write frequency file
            with open(kmerfr, "w") as outfile:

                for kmer, count in counts.most_common():
                    outfile.write(f"{kmer} : {count}\n")


            log("K-mers Frequency.txt is ready.")

            if ambiguous_bases > 0:
                log(
                    f"Warning: {ambiguous_bases} ambiguous bases detected. "
                    f"{skipped_kmers} k-mers were skipped."
                )

            hashUse(kmerfr, 0)

def runKmers(inputseq, k, generate_frequency, generate_all_kmers):
    
    if not generate_frequency and not generate_all_kmers:
        print("Can't perform k-mers algorithm. Please check the given parameters (Both are False).")
        log("Both runKmers parameters are False, terminating.")
        return
    
    print(">> Running K-mers for: ", inputseq.name, "\n")
    log("--Kmers algorith started for", inputseq.name, "--") 
    
    #get bin
    BASE_DIR = Path(__file__).resolve().parents[1] 

    #Get Kmers folder
    KmersDir = BASE_DIR / "Bin" / "Kmers" 
    
    #Get Kmers folder name
    KmersDirName = KmersDir / inputseq.stem
    
    
    '''
    Depending on why the algorithm is called, the parameters `generate_frequency`
    and `generate_all_kmers` determine which operations should be performed.
    
    During the initial call (from `main`), these parameters contain the user's
    requested actions. However, the algorithm checks which output files already
    exist and automatically adjusts the boolean values.
    
    This prevents unnecessary recalculation and allows each file for each organism
    to be processed only once, significantly improving the overall performance.
    '''
    
    if generate_frequency and generate_all_kmers: #in case both parameters are true 
        generate_frequency, generate_all_kmers = checkBothFiles(KmersDirName)  
    else:       #in case one of two parameters is True
        if checkExisting(KmersDirName,generate_frequency, generate_all_kmers):
            return # files already exist and are valid. No need to run k-mers again 

    #files doesn't exist or they are corrupted -> k-mers will run again for the given genome. 
    genomeDir = KmersDir / inputseq.stem
    genomeDir.mkdir(parents=True, exist_ok=True)
    
    if generate_frequency: 
        print("Counting genome's k-mers frequency.")
        log("Running generate_frequency")
        kmerfr = genomeDir / ("freqKmer_"+inputseq.stem+".txt")
        countkmerfrequency(inputseq,k,kmerfr)
    
    if generate_all_kmers:
        print("Creating all-kmers file.")
        log("Running generate_all_kmers")
        everykmer = genomeDir / ("allKmer_"+inputseq.stem+".txt")
        getKmers(inputseq,k,everykmer)
    
    print("\nK-mers finished.")
    log("--Kmers algorithm completed--\n")

def checkExisting(KmersDirName, generate_frequency, generate_all_kmers):
    if KmersDirName.exists() and KmersDirName.is_dir():                   #ckeck if the dir exist
        if checkFiles(KmersDirName, generate_frequency, generate_all_kmers):  #ckeck if the corresponding files (dipending the boolean parameters) do exist
            #if yes: 
            log("Κ-mers files for", KmersDirName.name, "already exist")
            log("--Kmers algorithm terminated.--\n")
            print("Κ-mers files for", KmersDirName.name, "already exist")
            print("--Kmers algorithm terminated.--\n")
            return True
        else:
            log("Couldn't locate the k-mers files for this genome. > Creating new ones.")
            print("Couldn't locate the k-mers files for this genome. > Creating new ones.")
            return False
    else:
        return False
    
    
def checkBothFiles(path):
    
    '''
    When both parameters are set to True, the program checks whether the required files already exist.
    If one or both files are found, the corresponding boolean values are updated so that only the necessary
    parts of the algorithm are executed.
    '''
    
    needFreqKmer = False
    needAllKmer = False
    
    if (path.exists() and path.is_dir()):  #ckeck if the dir exist
    
        if not any(path.iterdir()): # check if dir is empty. If yes, no need to continue, run the algorith anyway for both parameters. 
            log("Couldn't locate the k-mers files for this genome. > Creating new ones.")
            print("Couldn't locate the k-mers files for this genome. > Creating new ones.")
            return True, True
        
        filename1 = "freqKmer_"+path.stem+".txt"
        filepath1 = path / filename1
    
        filename2 = "allKmer_"+path.stem+".txt"
        filepath2 = path / filename2
        
        if not (filepath1.exists() and filepath1.is_file() and hashUse(filepath1, 1)):
            needFreqKmer = True
            
        if not (filepath2.exists() and filepath2.is_file() and hashUse(filepath2, 1)):
            needAllKmer = True
    else: 
        log("Couldn't locate the k-mers files for this genome. > Creating new ones.")
        print("Couldn't locate the k-mers files for this genome. > Creating new ones.")
        needFreqKmer = needAllKmer = True
    
    if not needFreqKmer and not needAllKmer:
        log("Κ-mers files for", path.stem, "already exist")
        print("Κ-mers files for", path.stem, "already exist")
        
    return needFreqKmer, needAllKmer
    
def checkFiles(path, generate_frequency, generate_all_kmers):
    
    filename1 = "freqKmer_"+path.stem+".txt"
    filepath1 = path / filename1

    filename2 = "allKmer_"+path.stem+".txt"
    filepath2 = path / filename2

    if not any(path.iterdir()): # check if dir is empty. If yes, no need to continue, return false. 
        return False
    
    if generate_frequency and not generate_all_kmers:
        if filepath1.exists() and filepath1.is_file() and hashUse(filepath1, 1): # check if files exist. If yes (they do exist), check if they are corrupted (hashUse).
            return True                                                          #if they exist and not corupted return True                                                      
        
    if generate_all_kmers and not generate_frequency:
        if filepath2.exists() and filepath2.is_file() and hashUse(filepath2, 1):
            return True
        
    return False                # Otherwise return false (meaning: files didn't exist or they are corrupted)
    
    
    
    
    
    
    
