import os
from collections import Counter
from pathlib import Path
from log import log
from hash_utils import hashUse
import json
from datetime import datetime
from progressBar import ProgressBar
from Stats import genomeTotalLines

version = "4.0"

'''
About this version: 
In this version, metadata files are generated for the output files.
Metadata have many practical uses, but in this implementation they are mainly
used to compare existing files with the ones that are about to be generated.
Also, progress bars have been added in order to track the analysis progress for larger genomes. 
------------

About k-mers algorithm:
The k-mer algorithm extracts all possible subsequences of length k from a DNA sequence.
Each k-mer represents a short fragment of the genome and can be used for sequence
comparison, similarity analysis, and genome characterization.
'''

def getKmers(inputseq,k,everykmer): 
    log("-GetKmers started-")
    
    progress = ProgressBar(
    total = genomeTotalLines(inputseq),
    description="Exporting all k-mers if file: ")
    
    kmercounter=0
    if ((os.path.getsize(inputseq) != 0) and (os.path.getsize(inputseq) >= k)): #If file is not empty, and if it is greater than k, proceed,
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
                prev = fulline[-(k - 1):] if k > 1 else "" 
                progress.update()
                
        hashUse(everykmer, 0)     #create sha256 file to secure authenticity
    else:
        log("Error while trying to read the file.")
    log(kmercounter, "k-mers of size", k , "have been found.")
    log("-GetKmers finished-")

def countKmerFrequency(inputseq,k,kmerfr):
    
            """
            Reads a genome file and creates a k-mer frequency file.

            The function:
            1) Reads the genome directly without creating an intermediate k-mer file (unless this action is needed -> getKmers function).
            2) Validates nucleotide symbols.
            3) Counts k-mer occurrences.
            4) Creates an output file containing k-mers sorted by frequency.
            """
            
            log("-countKmerFrequency started-")
            valid_bases = set("ATCGN")
            counts = Counter()   # Creates an empty Counter object, which is a specialized dictionary
                                 # from the collections module used for counting occurrences of elements.
                               
            ambiguous_bases = 0  # Counts the total number of ambiguous nucleotide symbols (e.g. 'N') found in the input genome.                   
            skipped_kmers = 0    # Counts the number of generated k-mers that were skipped because
                                 # they contain one or more ambiguous nucleotide symbols.

                                              
            progress = ProgressBar(
            total = genomeTotalLines(inputseq),
            description="Exporting k-mer frequencies: "
            )
                                 
            # Check if file exists and is large enough
            if not os.path.exists(inputseq):
                raise FileNotFoundError(f"Genome file not found: {inputseq}")        # Raises an exception and stops execution if the input genome file does not exist.
                log("Genome file not found: {inputseq}")
            if (os.path.getsize(inputseq) < k) or (os.path.getsize(inputseq) == 0):
                raise ValueError("Genome file is empty or smaller than given k.")    # Raises an exception and stops execution if the input genome file is empty.
                log("Genome file is empty or smaller than given k.")
            
            prev = "" # Stores the last (k-1) bases from the previous sequence line
            
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
            
                    prev = full_sequence[-(k - 1):] if k > 1 else "" # Store the last (k-1) bases for the next iteration
                    progress.update()

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
            log("-countKmerFrequency finished-")
            return ambiguous_bases
        

def runKmers(inputseq, k, generate_frequency, generate_all_kmers):
    '''
    Depending on why the algorithm is called, the parameters `generate_frequency`
    and `generate_all_kmers` determine which operations should be performed.
    
    During the initial call (from `main`), these parameters contain the user's
    requested actions. However, the algorithm checks which output files already
    exist and automatically adjusts the boolean values.
    
    This prevents unnecessary recalculation and allows each file for each organism
    to be processed only once, significantly improving the overall performance.
    '''
    
    log("Kmers version: ", version)
    
    if not generate_frequency and not generate_all_kmers:
        print("Can't perform k-mers algorithm. Please check the given parameters (Both are False).")
        log("Both runKmers parameters are False, terminating.")
        return
    

    print("-Running K-mers for:", inputseq.name,"-\n")
    log("-runKmers started for", inputseq.name, "-") 
    
    #get bin
    BASE_DIR = Path(__file__).resolve().parents[1] 

    #Get Kmers folder path
    KmersDir = BASE_DIR / "Bin" / "Kmers" 
    
    #Get Kmers folder path with dir name 
    KmersDirName = KmersDir / inputseq.stem #.stem cuts the extension
    
    '''
    The below if-else statement validates requested k-mer outputs through metadata comparison.
    Regenerates only missing or incompatible files, preventing redundant k-mer calculations
    when valid results already exist.
    '''
    
    if generate_frequency and generate_all_kmers:    #in case both parameters are true 
    
        #read metadata for both parameters 
        #if k in metadata files and given k arg agree, continue, else run kmers for new boolean values
                                                                       
        metaData = readMetadata(KmersDirName, "frequency")  #first read metadata       
        
        if metaData["output_type"] == "idle" or int(metaData["k"]) != k: #if metadata don't exist/is corrupted (idle mode) or exist but for different k than the given k
            runKmers(inputseq, k, True, False)                           #run the algorith with the right parameters

        metaData = readMetadata(KmersDirName, "allKmers")                
        
        if metaData["output_type"] == "idle" or int(metaData["k"]) != k: #same check as above but for the second file 
            runKmers(inputseq, k, False, True)                           #(remember, we are in the case (outer if) where the requested use includes both parameters)
        
        generate_frequency, generate_all_kmers = checkBothFiles(KmersDirName)  #if metadata pass the checks, then check also the generated files (that should exist at this point)
        
        if not generate_frequency and not generate_all_kmers:                           #if bothe boolean values turns out to be False then,
            return                                                                      #the corresponding files already exist and are valid. No need to run k-mers again 
    else: #one of the two parameters is True and the other is False                                                                               #in case one of two parameters is True
        #first check the metadata for each case
        
        skipCheck = False #this is a boolean flag variable that determines if there's need to run the check
        
        if generate_frequency and not generate_all_kmers: #True, False
        
            metaData = readMetadata(KmersDirName, "frequency")  
            if metaData["output_type"] == "idle" or int(metaData["k"]) != k: #same logic as above
                generate_frequency = skipCheck = True
        else: #opposite case
        
            metaData = readMetadata(KmersDirName, "allKmers")
            if metaData["output_type"] == "idle" or int(metaData["k"]) != k:
                generate_all_kmers = skipCheck = True
       
        '''
        Logic behind "skipCheck":
        skipCheck avoids unnecessary file validation when metadata indicates that existing files
        cannot be reused (idle state or different k value).
        
        In this case, the corresponding generation flag is set to True so the required k-mer
        calculation will run later and produce new valid output files. Therefore, calling
        checkExisting() is unnecessary because the existing files will be regenerated.
        '''

        if checkExisting(KmersDirName,generate_frequency, generate_all_kmers) and not skipCheck: 
            return  # files already exist and are valid. No need to run k-mers again (return)
        

    #if files doesn't exist or they are corrupted -> k-mers will run again for the given genome.
    
    genomeDir = KmersDir / inputseq.stem
    genomeDir.mkdir(parents=True, exist_ok=True)
    
    if generate_frequency: 
        print("Counting genome's k-mers frequency.")
        log("Running generate_frequency")
        kmerfr = genomeDir / ("freqKmer_"+inputseq.stem+".txt")
        countKmerFrequency(inputseq,k,kmerfr)
        generateMetadata(KmersDirName, k, "frequency")
    
    if generate_all_kmers:
        print("Creating all-kmers file.")
        log("Running generate_all_kmers")
        everykmer = genomeDir / ("allKmer_"+inputseq.stem+".txt")
        getKmers(inputseq,k,everykmer)
        generateMetadata(KmersDirName, k, "allKmers")

    print("\nK-mers finished.")
    log("-runKmers finished-")

def checkExisting(KmersDirName, generate_frequency, generate_all_kmers):
    if KmersDirName.exists() and KmersDirName.is_dir():                       #ckeck if the dir exist
        if checkFiles(KmersDirName, generate_frequency, generate_all_kmers):  #ckeck if the corresponding files (dipending the boolean parameters) exist
            #if yes: 
            log("Κ-mers files for", KmersDirName.name, "already exist")
            print("Κ-mers files for", KmersDirName.name, "already exist")
            return True
        else:
            log("Missing k-mers files for this genome. > Creating new ones.")
            print("Missing k-mers files for this genome. > Creating new ones.")
            return False
    else:
        return False
    
    
def checkBothFiles(path):
    
    '''
    When both parameters are set to True, the program checks whether the required files already exist.
    If one or both files are found, the corresponding boolean values are updated so that only the necessary
    parts of the algorithm are executed.
    '''
    needFreqKmer = needAllKmer = False
    
    if (path.exists() and path.is_dir()):  #ckeck if the dir exist
        if not any(path.iterdir()): # check if dir is empty. 
            log("Missing k-mers files for this genome. > Creating new ones.")
            print("Missing k-mers files for this genome. > Creating new ones.")
            return True, True #If yes, no need to continue, run the algorith anyway for both parameters. 
        
        filename1 = "freqKmer_"+path.stem+".txt"
        filepath1 = path / filename1
    
        filename2 = "allKmer_"+path.stem+".txt"
        filepath2 = path / filename2
        
        if not (filepath1.exists() and filepath1.is_file() and hashUse(filepath1, 1)):
            needFreqKmer = True
            
        if not (filepath2.exists() and filepath2.is_file() and hashUse(filepath2, 1)):
            needAllKmer = True
    else: 
        log("Missing k-mers files for this genome. > Creating new ones.")
        print("Missing k-mers files for this genome. > Creating new ones.")
        needFreqKmer = needAllKmer = True
    
    if not needFreqKmer and not needAllKmer:
        log("Κ-mers files for", path.stem, "already exist")
        #print("Κ-mers files for", path.stem, "already exist")
        
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
    

def generateMetadata(source_genome, k, output_type):
    
    """
    Creates a metadata JSON file for generated k-mer files.

    Metadata stores information about the generated file and
    the parameters used during its creation.
    """
    
    metaDataFile = source_genome.name+"_"+output_type+"_metaData.json"
    metaDataFilePath = source_genome / metaDataFile
    
    metadata = {
        "Version": version,
        "source_genome": source_genome.name,
        "output_file": metaDataFile,
        "output_type": output_type,
        "k": k,
        "created_at": datetime.now().isoformat()
    }
    
    if output_type != "idle": 
        with open(metaDataFilePath, "w") as outfile:
            json.dump(metadata, outfile, indent=4) # The indent parameter determines the number of spaces used for indentation
                                                   # when formatting the JSON files
        hashUse(metaDataFilePath, 0)               #create sha256 file to secure authenticity
        
    return metadata
    
def readMetadata(source_genome, output_type):
    #if the corresponding metadata doesn't exist or sha fails return false
    """
    Reads a metadata JSON file and returns its contents as a dictionary.
    """
    
    metaDataFile = source_genome.name+"_"+output_type+"_metaData.json"
    metaDataFilePath = source_genome / metaDataFile 
    
    
    #check if file exist and if it's valid
    if metaDataFilePath.exists() and metaDataFilePath.is_file() and hashUse(metaDataFilePath, 1):
        #metadata file exist and it's not corrupted
        with open(metaDataFilePath, "r") as infile:
            metadata = json.load(infile)
            return metadata
    else:
        return generateMetadata(source_genome, 0, "idle")
        
    
    

   
    
