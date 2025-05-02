#!/usr/bin/python
import sys
def get_ids(idlist):
    # Read the file containing the list of IDs (one per line) and return them as a list
    f = open(idlist)
    return f.read().strip().split('\n')

def get_seq(pidlist, seqfile):
    f = open(seqfile)
    s = 0
    for line in f:
        # The FASTA format is: a line starting with '>' (the header), followed by the sequence on the next line(s)
        # We define a state variable `s` to control whether to print the current sequence or not
        if line.startswith('>'):
            pid = line.split('|')[1]  # Extract the accession number from the second field of the header
            s = 0  # Default state: skip the sequence
            if pid in pidlist:
                s = 1  # If the accession number is in the list, enable printing
        if s == 1:
            print(line.strip())

# NOTE:
# The input FASTA file must follow the UniProt format:
# >sp|P12345|... (ID in 2nd field)
# Before using this script, ensure that your ID list corresponds to the accession number field



if __name__ == '__main__':
    idlist = sys.argv[1]  # File containing the list of IDs to keep
    seqfile = sys.argv[2]  # FASTA file from which to extract the matching sequences
    pidlist = get_ids(idlist)
    get_seq(pidlist, seqfile)