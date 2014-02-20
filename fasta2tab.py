#!/usr/bin/env python
import random
import string
from Bio import SeqIO

# the character set for sequence UUIDs
charset = string.letters + string.digits
# the number of characters to use for sequence UUIDs
n_char = 8

def fasta2tab_file(fname):
    '''
    This function will process a single clip fasta file into a single tab file
    appropriately named.
    '''
def read_fasta_file(fname_fasta):
    '''
    This function will read a clip fasta file, and return a list of some kind
    of Bio.SeqRecord that will contain all of the metadata about a sequence
    and the sequence itself.
    Likely, it will use attribute names like:
    SEQUENCE_ID
    SEQUENCE_GAP
    GERMLINE_GAP_DMASK
    CLONE
    '''
    f = open(fname_fasta,'r')
    lst_seq_record = [record for record in SeqIO.parse(f,'fasta')]
    n_id = len(lst_seq_record)
    lst_uuid = get_n_uuid(n_id)
    
    return lst_seq_record

def write_tab_file(fname_fasta,lst_seq_record):
    '''
    This function will take a list of Bio.SeqRecord instances, extract the
    needed sequence and metadata from each, and write them to an appropriately
    formatted tab file.
    '''
    with open(fname_fasta,'w') as f:
        tab_writer = csv.writer(f, delimiter='\t')
        # write header
        for seq_record in lst_seq_record:
            # write row
#            tab_writer.writerow([
#                seq_record.SEQUENCE_ID,
#                seq_record.SEQUENCE_GAP,
#                seq_record.GERMLINE_GAP_DMASK,
#                seq_record.CLONE,
#                ])
            pass

def get_n_uuid(n):
    '''
    This function will generate as many UUID as requested, confirming that
    they are unique.
    It is neither elegant nor efficient.
    '''
    # generates list of UUID prior to checking for duplicates
    lst_uuid = [get_uuid() for iI in xrange(n)]
    # WARNING: while loop, remove in future versions
    # check to see that all UUID in lst_uuid are unique
    while len(lst_uuid) > len(set(lst_uuid)):
        lst_uuid = [get_uuid() for iI in xrange(n)]
    return lst_uuid
def get_uuid():
    '''
    This function returns a single UUID using module variables charset,n_char.
    '''
    return ''.join(random.choice(charset) for iI in xrange(n_char))
if __name__ == '__main__':
    import sys
    lst_fasta_files = sys.argv[1:]
    for fname in lst_fasta_files:
        fasta2tab_file(fname)

