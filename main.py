# Final project


'''
Process Frow:
De Bruijn Garph -> Eulerian Cycle/Path -> Alignment
     (1)        ->         (2)         ->    (3)

(1) FASTQ file of Reads -> Graph PDF, edges.txt
'''

import sys
import os
import numpy as np
import de_bruijn as db
import alignment as al
import eulerianCycle as ec
import eulerianPath as ep

# ==============================================================================================================
# Retreives header/sequence pair data from the fna file
# Retruns seq_data: list of sequences, header_data: List of corresponding header info
def get_data(filename):

    seq_data = []
    header_data = []
    mime = filename.split('.').pop()                        # Get file MIME type

    if(mime == 'fna'):                                      # Input file is a fna (data every two lines)
        print("Reading FNA file: ", filename)

        with open(filename) as f:
            for head, seq in zip(f,f):                      # Get header and sequence info
                head, seq = head.strip(), seq.strip()       # Strip newline characters
                header_data.append(head[1:])                # Add to header list
                seq_data.append(seq)                        # Add to sequece list

    elif(mime == 'fastq'):                                  # File is a fastq (data every four lines)
        print("Reading FASTQ file: ", filename)

        with open(filename) as f:
            for head, seq, p, score in zip(f,f,f,f):        # Get four line at a time (Header, sequence, plus thingy, score)
                head, seq = head.strip(), seq.strip()       # Strip header and sequece data of newline chars
                header_data.append(head[1:])                # Add to header list
                seq_data.append(seq)                        # Add to sequence list

    elif(mime == 'txt'):                                    # File is a Text
        print("Reading text file ", filename)

        with open(filename) as f:
            for seq in f:                                   # Read each line
                seq = seq.strip()                           # Strip data of new line characters
                seq_data.append(seq)                        # Append data to list

            header_data.append('Assembled data')            # No headers should be in the file so add this one
    else:
        print("ERROR: Invalid File Type!")
        print("Only fna, fastq, or txt types! The entered file type is ", mime)
        return

    return (seq_data, header_data)

# ==============================================================================================================
# Writes data to a text file
def make_txt(data, filename='./output/temp/output.txt'):

    with open(filename, 'w', newline='') as file:
        for x in data:
            file.write(str(x) + '\n')

# ==============================================================================================================
# Loops thru a k-mer range, builds a set of kmers and edges
def loop_kmer(data, kstart, kend, lstart=0, lend=1):
    db_data = db.De_bruijn(data[0], data[1])                    # Create de Bruijn graph

    for i in range(kstart, kend+1):                             # Loop thru k-mer range
        db_data.de_bruijn_graph(k=i)
        db_data.make_docs(True, True, True, str(i))

# ==============================================================================================================
# Clears out old files
def rmove():
    dir = './output/align/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in align directory

    dir = './output/eulerian/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in eulerian directory

    dir = './output/graph/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in graph directory

    dir = './output/temp/'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))     # Remove files in temp directory

# ==============================================================================================================
# Finds eulerian cycle and path, writes the data to a text file
def eulerian_string(filename, kmer):

    text = []
    with open(filename, 'r') as f:                              # Open directed graph file
        for line in f:
            text.append(line.strip())                           # Add the edges and strip the newline characters

    cycle = ec.EulerianCycle(text)                              # Find eulerian cycle
    path = ep.EulerianPath(text)                                # Find eulerian cycle
    cycle_file = 'output/eulerian/eulerianCycle_' + kmer + '.txt'  # File that will be written to
    path_file = 'output/eulerian/eulerianPath_' + kmer + '.txt'    # File that will be written to

    with open(cycle_file, 'w') as f:                            # Write eulerian cycle to file
        f.write(''.join(cycle))

    with open(path_file, 'w') as f:                             # Write eulerian path to file
        f.write(''.join(path))

# ==============================================================================================================
# Aligns assembled contig with reference genome, creates text file of alignment, and plots comparison
def align_contig():
    dir = './output/align/'                                         # Directory to beused
    ref_seg = get_data('./input/sars_spike_protein_assembled.fna')  # Get data for reference sequence (Spike protien)

    for file in os.listdir(dir):                                    # Iterate thru files in align folder
        f = os.path.join(dir, file)
        print(f)
        data = get_data(f)                                          # Get in data from files for alignment

        a = al.alignment(ref_seg[0], data[0])                       # Initialize alignment
        a.get_alignment()                                           # Retrieve alignment data
        #a.plot_comp(filename)                                      # Create a comparison plot of the two sequences

# ==============================================================================================================
def main():

    # ----------------------------------------------------------------------------------------------------------
    # User enters prefered file
    if(len(sys.argv) == 2):
        rmove()
        inpt = sys.argv[1]                                          # Extracting fastq file name from arguments
        file = "./input/sars_spike_protein_reads.fastq"

        if(inpt == '-l'):                                           # User enters one kmer value and read range
            k = int(input("Enter k-mer value: "))
            lstart = int(input("Enter index for first sequence: "))
            lend = int(input("End index for the last sequence: "))

            data = get_data(file)                                   # Processing genome data from fna file
            db_graph = db.De_bruijn(data[0], data[1])               # Create de Bruijn graph
            db_graph.de_bruijn_graph(start=lstart, end=lend, k=k)
            db_graph.make_docs(True,True,True,str(k))

        elif(inpt == '-k'):                                           # User enters kmer range
            kstart = int(input("Enter starting k-mer value: "))
            kend = int(input("Enter ending k-mer value: "))

            data = get_data(file)                                   # Processing genome data from fna file
            loop_kmer(data, kstart, kend)

        elif(inpt == '-kl'):                                        # User enters kmer range and read range
            kstart = int(input("Enter starting k-mer value: "))
            kend = int(input("Enter ending k-mer value: "))
            lstart = int(input("Enter index for first sequence: "))
            lend = int(input("End index for the last sequence: "))

            data = get_data(file)                                   # Processing genome data from fna file
            loop_kmer(data, kstart, kend, lstart, lend)

        elif(inpt == '-a'):                                         # Align contig with reference genome (Spike protein)
            data = get_data(file)
            align_contig(data)

        else:
            print("ERROR: Invalid Input!")
            print("Use:\n-l\tsegmenting reads list\n-k\tk-mer range\n-kl\tk-mer range and segmenting reads list\n-a\taligning contig")
            return
    # ----------------------------------------------------------------------------------------------------------
    # Use default file
    elif(len(sys.argv) == 1):
        rmove()

        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                       # Processing genome data from fna file
        k = 10

        db_graph = db.De_bruijn(data[0], data[1])                       # Create de Bruijn graph
        db_graph.de_bruijn_graph(start=0, end=2, k=k)
        db_graph.make_docs(True,True,True, str(k))

    # ----------------------------------------------------------------------------------------------------------
    # Use default file and run kmer loop
    # python .\main start stop
    '''
    elif(len(sys.argv) == 3):
        rmove()

        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                       # Processing genome data from fna file

    '''
if __name__ == "__main__":
    main()
