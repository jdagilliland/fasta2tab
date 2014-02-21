#!/usr/bin/env python
import random
import string
import csv

from Bio import SeqIO

# the character set for sequence UUIDs
charset = string.letters + string.digits
# the number of characters to use for sequence UUIDs
n_char = 8

def prune_germline_records(lst_seq_record):
    return [record for record in lst_seq_record if
        record.description[0] !='>']

def parse_clip_attr(seq_record):
    '''
    This function accepts a Bio.SeqRecord, and parses the
    'description' into a dictionary of supplementary 'clip' attributes, with
    which it then updates the attributes of the SeqRecord.
    '''
    # splits description by '|' into key-value pair strings
    lst_str_keyval = [keyval for keyval in seq_record.description.split('|')]
    # splits each key value pair string by ':', and creates dictionary
    try:
        dict_clip_attr = dict([keyval.split(':',1) for keyval in lst_str_keyval
            if len(keyval)>0])
    except:
        print(lst_str_keyval)
        raise
    # updates seq_record with those dictionary items
    seq_record.__dict__.update(dict_clip_attr)
    return seq_record

def append_uuid(lst_seq_record):
    '''
    This function accepts a list of Bio.SeqRecord, computes a unique
    identifier for each, and appends it to the 'SEQUENCE_ID' string of each
    SeqRecord.
    '''
    n_id = len(lst_seq_record)
    lst_uuid = get_n_uuid(n_id)
    for iI, record in enumerate(lst_seq_record):
        try:
            record.SEQUENCE_ID += lst_uuid[iI]
        except:
            print(record)
            raise
    
def tabname(fastaname):
    '''
    This function generates a suitable name for the tab file based on the name
    of the fasta file.
    '''
    return fastaname.rpartition('.')[0] + '.tab'

def read_fasta_file(fname_fasta):
    '''
    This function will read a clip fasta file, and return a list of some kind
    of Bio.SeqRecord that will contain all of the metadata about a sequence
    and the sequence itself.
    Likely, it will use attribute names like:
    SEQUENCE_ID
    SEQUENCE_GAP
    GERMLINE_GAP_DMASK
    CLONE_ID
    '''
    f = open(fname_fasta,'r')
    lst_seq_record = [record for record in SeqIO.parse(f,'fasta')]
    return lst_seq_record

tpl_cols = (
    'SEQUENCE_ID',
    'CLONE_ID',
    )
def write_tab_file(fname_tab,lst_seq_record):
    '''
    This function will take a list of Bio.SeqRecord instances, extract the
    needed sequence and metadata from each, and write them to an appropriately
    formatted tab file.
    '''
    with open(fname_tab,'w') as f:
        tab_writer = csv.writer(f, delimiter='\t')
        # write header
        tab_writer.writerow(tpl_cols)
        for seq_record in lst_seq_record:
            # write row
            tab_writer.writerow([seq_record.__dict__[field] for field
                in tpl_cols])
            pass

def get_n_uuid(n):
    '''
    This function will generate as many UUID as requested, confirming that
    they are unique.
    '''
    # generates list of UUID prior to checking for duplicates
    lst_uuid = [get_uuid() for iI in xrange(n)]
    # WARNING: while loop, remove in future versions
    # check to see that all UUID in lst_uuid are unique
    while True:
        lst_uuid = set(lst_uuid)
        n_unique = len(lst_uuid)
        if n_unique == n:
            break
        lst_uuid.update([get_uuid() for iI in xrange(n-n_unique)])
    return list(lst_uuid)

def get_uuid():
    '''
    This function returns a single UUID using module variables charset,n_char.
    '''
    return ''.join(random.choice(charset) for iI in xrange(n_char))

def fasta2tab_file(fname):
    '''
    This function will process a single clip fasta file into a single tab file
    appropriately named.
    '''
    #generate an initial list of SeqRecord
    lst_seq_record = read_fasta_file(fname)
    #prune out entries that list germlines
    lst_seq_record = prune_germline_records(lst_seq_record)
    #parse supplementary attributes
    for record in lst_seq_record:
        parse_clip_attr(record)
    #append uuid portion to SEQUENCE_ID
    append_uuid(lst_seq_record)
    #find an appropriate name for the tab file
    tabfname = tabname(fname)
    #write the final list of SeqRecord to the tab file
    write_tab_file(tabfname, lst_seq_record)
    return None

if __name__ == '__main__':
    import sys
    lst_fasta_files = sys.argv[1:]
    for fname in lst_fasta_files:
        fasta2tab_file(fname)

