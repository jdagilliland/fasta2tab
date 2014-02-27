#!/usr/bin/env python
import random
import string
import csv

from Bio import SeqIO

# the character set for sequence UUIDs
charset = string.letters + string.digits
# the number of characters to use for sequence UUIDs
n_char = 8
# these are the columns to output in the tab file
tpl_cols = (
    'SEQUENCE_ID',
    'SEQUENCE',
    'FUNCTIONAL',
    'IN_FRAME',
    'STOP',
    'MUTATED_INVARIANT',
    'INDELS',
    'V_MATCH',
    'V_LENGTH',
    'J_MATCH',
    'J_LENGTH',
    'V_CALL',
    'D_CALL',
    'J_CALL',
    'SEQUENCE_GAP',
    'V_GAP_LENGTH',
    'N1_LENGTH',
    'D_5_TRIM',
    'D_3_TRIM',
    'N2_LENGTH',
    'J_5_TRIM',
    'J_GAP_LENGTH',
    'JUNCTION_GAP_LENGTH',
    'JUNCTION',
    'PRIMER',
    'CONSCOUNT',
    'DUPCOUNT',
    'CLONE',
    'GERMLINE_GAP_D_MASK',
    )

class ClipRecord:
    mask_char = 'n'
    def __init__(self, seq_record):
        self.seq_record = seq_record
        self.description = self.seq_record.description
        self.seq = self.seq_record.seq
        self.SEQUENCE = str(self.seq_record.seq).lower()
        # splits description by '|' into key-value pair strings
        lst_str_keyval = [keyval for keyval in self.description.split('|')]
        # splits each key value pair string by ':', and creates dictionary
        try:
            dict_clip_attr = dict([keyval.split(':',1) for keyval in lst_str_keyval
                if len(keyval)>0])
        except:
            print(lst_str_keyval)
            raise
        for key, value in dict_clip_attr.items():
            setattr(self, key, value)
        if hasattr(self, 'CLONE_ID'):
            self.CLONE = self.CLONE_ID
        if hasattr(self, 'SEQUENCE'):
            self.SEQUENCE_GAP = self.SEQUENCE
        return None
    def d_mask(self):
        '''
        This method takes the sequence data, uses the fields indicating
        V_LENGTH and J_LENGTH, and masks what is in between (the D region)
        with class attribute: mask_char, outputting to field
        GERMLINE_GAP_D_MASK.
        '''
        # check to make sure that the ClipRecord instance has all the necessary fields
        #if not has_attr(self,...)
        v_length = int(self.V_LENGTH)
        j_length = int(self.J_LENGTH)
        length = len(self.SEQUENCE)
        d_length = length - v_length - j_length
        lst_seq = list(self.SEQUENCE)
        
        if len(lst_seq) != length:
            print('Length prior to masking')
            print(length)
            print('Length after masking')
            print(len(lst_seq))
            print(v_length)
            print(j_length)
            print(d_length)
#            raise ValueError('''The sequence has changed in length.
#                You must mask a different way.''')
        self.GERMLINE_GAP_D_MASK = ''.join(lst_seq)
        pass
    def d_no_mask(self):
        '''
        This method takes the sequence data and copies it into
        GERMLINE_GAP_D_MASK without actually masking the D region.
        '''
        self.GERMLINE_GAP_D_MASK = self.SEQUENCE
        pass
    pass
def prune_germline_records(lst_seq_record):
    lst_germline = [record for record in lst_seq_record if
        record.description[0] =='>']
    lst_non_germline = [record for record in lst_seq_record if
        record.description[0] !='>']
    return lst_germline, lst_non_germline

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
    of object that will contain all of the metadata about a sequence
    and the sequence itself.
    Likely, it will use attribute names like:
    SEQUENCE_ID
    SEQUENCE_GAP
    GERMLINE_GAP_DMASK
    CLONE_ID
    '''
    f = open(fname_fasta,'r')
    lst_seq_record = [ClipRecord(record) for record in SeqIO.parse(f,'fasta')]
    
    return lst_seq_record

def write_tab_file(fname_tab,lst_seq_record):
    '''
    This function will take a list of Bio.SeqRecord instances, extract the
    needed sequence and metadata from each, and write them to an appropriately
    formatted tab file.
    '''
    with open(fname_tab,'w') as f:
        tab_writer = csv.DictWriter(f, tpl_cols, extrasaction='ignore', delimiter='\t')
        # write header
        tab_writer.writeheader()
        for seq_record in lst_seq_record:
            tab_writer.writerow(seq_record.__dict__)
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
    lst_germline, lst_clip_seq = prune_germline_records(lst_seq_record)
    #parse supplementary attributes
    #mask d region
    for clip_seq in lst_clip_seq:
        clip_seq.d_no_mask()
    #append uuid portion to SEQUENCE_ID
    append_uuid(lst_clip_seq)
    #find an appropriate name for the tab file
    tabfname = tabname(fname)
    #write the final list of SeqRecord to the tab file
    write_tab_file(tabfname, lst_clip_seq)
    return None

if __name__ == '__main__':
    import sys
    lst_fasta_files = sys.argv[1:]
    for fname in lst_fasta_files:
        fasta2tab_file(fname)

