# Final project

import sys
import numpy as np
import de_bruijn as db

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

    else:
        print("ERROR: Invalid File Type!")
        print("Only fna or fastq types, entered file type is ", mime)
        return

    return (seq_data, header_data)

# ==============================================================================================================
# Writes data to a text file
def make_txt(data, filename='./output/output.txt'):

    with open(filename, 'w', newline='') as file:
        for x in data:
            file.write(str(x) + '\n')


def create_directed_graph(edges):
    with open("output/spike_protein_directed_graph.txt", "w") as f:
        added_nodes = set()                                 # Set of (nodes, destination) pairs that have already been iterated through

        for edge in edges:
            node,dest = edge
            if edge not in added_nodes:                     # Check if (node, destination) pair has already been writted
                f.write(node + ' -> ' + dest)

                for edge2 in edges:                         # Iterate through all edges to see if there are any other edges coming out of node
                    node2, dest2 = edge2
                    # If the nodes are the same AND (node, destination) pair is not the same as above AND this node has not already been created
                    if node == node2 and edge != edge2 and edge2 not in added_nodes:
                        f.write(',' + dest2)                # Write the additional destination for exising node
                        added_nodes.add(edge2)              # add (node, destination) pair to already added set

                f.write('\n')
    f.close()

# ==============================================================================================================
def main():

    # ----------------------------------------------------------------------------------------------------------
    # User enters prefered file
    if(len(sys.argv) == 2):
        fna_file = sys.argv[1]                                                          # Extracting fastq file name from arguments
        data = get_data(fna_file)                                                       # Processing genome data from fna file

        db_graph = db.De_bruijn(data[0], data[1])                                       # Create de Bruijn graph
        db_graph.matplot_graph(False,True)

    # ----------------------------------------------------------------------------------------------------------
    # Use default file
    elif(len(sys.argv) == 1):
        fna_file = "./input/sars_spike_protein_reads.fastq"
        data = get_data(fna_file)                                                       # Processing genome data from fna file

        db_graph = db.De_bruijn(data[0], data[1])                                       # Create de Bruijn graph
        #db_graph.matplot_graph(False,True)

        create_directed_graph(db_graph.edges)                                           # Create txt file containing directed graph of sars spike protein

if __name__ == "__main__":
    main()
