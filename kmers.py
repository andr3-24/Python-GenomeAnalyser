import os
from collections import Counter
from pathlib import Path
from log import log
from hash_utils import hashUse

'''
The k-mer algorithm extracts all possible subsequences of length k from a DNA sequence.
Each k-mer represents a short fragment of the genome and can be used for sequence
comparison, similarity analysis, and genome characterization.
'''

def checkFiles(path):
    filename1 = "allKmer_"+path.stem+".txt"
    filename2 = "freqKmer_"+path.stem+".txt"
    
    filepath1 = path / filename1
    filepath2 = path / filename2
    
    if not any(path.iterdir()):
        return False

    if filepath1.exists() and filepath1.is_file() and filepath2.exists() and filepath2.is_file():           # check if files exist
        if hashUse(filepath1, 1) and hashUse(filepath2, 1):                                                 # if yes (they do exist), check if they are corrupted.
            return True                                                                                     # if files pass the check, return true (meaning files exist and the are ready for use)
        
    return False                                                                                            # Otherwise return false (meaning: files didn't exist or they are corrupted)



def checkExisting(path, folder):
    existOrNot = path / folder
    if existOrNot.exists() and existOrNot.is_dir(): #ckeck if the dir exist
        if checkFiles(existOrNot): #if also the needed files already exist
            log("Κ-mers files for", folder, "already exist")
            log("--Kmers algorithm terminated.--\n")
            print("Κ-mers files for", folder, "already exist")
            print("--Kmers algorithm terminated.--\n")
        else:
            log("Couldn't locate the k-mers files for this genome. > Creating new ones.")
            print("Couldn't locate the k-mers files for this genome. > Creating new ones.")
            return False
    else:
        return False

def runKmers(inputseq, k):
    print(">> Running K-mers for: ", inputseq, "\n")
    log("--Kmers algorith started--")
    
    #get bin
    BASE_DIR = Path(__file__).resolve().parents[1] 

    #Get Kmers folder
    KmersDir = BASE_DIR / "Bin" / "Kmers" 
    
    #check if genome kmers folder already exists
    if checkExisting(KmersDir, inputseq.stem):
        return
    else:
        genomeDir = KmersDir / inputseq.stem
        genomeDir.mkdir(parents=True, exist_ok=True)
        
    
    everykmer = genomeDir / ("allKmer_"+inputseq.stem+".txt")
    kmerfr = genomeDir / ("freqKmer_"+inputseq.stem+".txt")

    filename = inputseq.name
    
    log("Analysing", Path(filename).stem, " genome.:")
    #print("Analysing ", Path(filename).stem, " genome.")
    
      
    def getKmers():
        kmercounter=0
        if ((os.path.getsize(inputseq) != 0) and (os.path.getsize(inputseq) >= k)): #If file is not empty (((os.path.getsize(inputseq) != 0)), and if it is bigger than k, proceed,
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
                
                
    def countkmerfrequency():
        distinct = set()
        
        with open(everykmer, "r") as inputfile, open(kmerfr, "w") as outfile:
            #this function uses the everykmer document and creates a new file such as (4)
            
            #1) get distinct values and place them in a data stracture
            #2) count the frequency of each kmer
            #3) create a new document containing the occurances of each distinct kmer, sorted in descending order. 
            
            #1
            for line in inputfile:
                distinct.add(line.rstrip("\n"))
                
            allkmersDistinct = list(distinct)
            log("There are "+ str(len(allkmersDistinct))+ " distinct k-mers.")
            #print("There are "+ str(len(allkmersDistinct))+ " distinct k-mers.")

            #2
            inputfile.seek(0) #Set inputfile index to 0
            counts = Counter(line.strip() for line in inputfile)
            
            #3
            for kmer, count in counts.most_common():
                outfile.write(f"{kmer} : {count}\n")
            
            log("Κ-mers Frequency.txt is ready.")
            #print(">> Κ-mers Frequency.txt is ready.")
            
        hashUse(kmerfr, 0)  #create sha256 file to secure authenticity 
        
    getKmers()
    countkmerfrequency()
    
    log("Analysing for ", Path(filename).stem , " is done. ")
    log("--Kmers algorithm completed--\n")
    
    print("K-mers finished.")
    #print("Analysing for ", Path(filename).stem , " is done. ")
            

