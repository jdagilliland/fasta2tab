#!/usr/bin/env python
import random
import string
import csv

from Bio import SeqIO

# the character set for sequence UUIDs
charset = string.letters + string.digits
# the number of characters to use for sequence UUIDs
n_char = 8
tpl_germ_cols = (
    'clone',
    'seq',
    )
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
class HeaderError(ValueError):
    pass
class ClipRecord:
    mask_char = 'n'
    field_sep = '|'
    colon = ':'
    def __init__(self, seq_record):
        self.seq_record = seq_record
        self.description = self.seq_record.description
        self.seq = self.seq_record.seq
        self.SEQUENCE = str(self.seq_record.seq).lower()
        # set fields from FASTA header row
        self._set_from_header()
        if hasattr(self, 'CLONE_ID'):
            if self.CLONE_ID[:9] == 'Germline:':
                self.CLONE = self.CLONE_ID[9:]
            else:
                self.CLONE = self.CLONE_ID
        if hasattr(self, 'SEQUENCE'):
            self.SEQUENCE_GAP = self.SEQUENCE
        if hasattr(self, 'INDELS'):
            if self.INDELS != 'T':
                self.INDELS = 'F'
        if not hasattr(self, 'INDELS'):
            self.INDELS = 'F'
        return None
    def _set_from_header(self):
        '''
        This method scrapes input fields from the FASTA sequence description
        in the format of '[|]key:value[|key:value]', and adds them all as
        instance attributes.
        This method is run by __init__, is private, because it should
        '''
        # keys_list will keep track of what keys are added from the header,
        # that way they can later be used to put them all in the TAB file.
        try:
            dict_clip_attr = parse_header_delim(
                self.description, self.field_sep, self.colon)
            # this lets other methods know which style the header was
            self.header_style = 0
        except HeaderError:
            # if parsing using the default scheme fails, try the one backup
            # scheme that I know might be used
            dict_clip_attr = parse_header_delim(
                self.description, '>', '|')
            # this lets other methods know which style the header was
            self.header_style = 1
        except:
            raise HeaderError
        self.keys_list = set()
        for key, value in dict_clip_attr.items():
            setattr(self, key, value)
            self.keys_list.update(key)
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
        
        self.GERMLINE_GAP_D_MASK = mask_d_region_length(
            self.SEQUENCE, v_length, j_length, mask_char=self.mask_char)
        pass
    def d_no_mask(self):
        '''
        This method takes the sequence data and copies it into
        GERMLINE_GAP_D_MASK without actually masking the D region.
        '''
        self.GERMLINE_GAP_D_MASK = self.SEQUENCE
        pass
    def germline_d_no_mask(self, dict_germline):
        '''
        This method populates the GERMLINE_GAP_D_MASK field with the
        appropriate sequence from the germline entries.
        '''
        self.GERMLINE_GAP_D_MASK = dict_germline[self.CLONE]
        pass
    def prep_germ(self):
        '''
        This method prepares a ClipRecord instance for entry into a germline
        TAB file.
        '''
        self.dict_germ_attr = dict()
        if hasattr(self,'Germline'):
            self.dict_germ_attr['clone'] = getattr(self,'Germline')
        elif hasattr(self,'>Germline'):
            self.dict_germ_attr['clone'] = getattr(self,'>Germline')
        else:
            raise HeaderError('''ClipRecord object has no suitable 'Germline' field.''')
        if hasattr(self,'SEQUENCE'):
            self.dict_germ_attr['seq'] = getattr(self,'SEQUENCE')
        else:
            raise HeaderError('''ClipRecord object has no suitable sequence.''')
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

def germtabname(fastaname):
    '''
    This function generates a suitable name for the germline tab file based
    on the name of the fasta file.
    '''
    return fastaname.rpartition('.')[0] + '_germ' + '.tab'

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
    lst_seq_record = [record for record in SeqIO.parse(f,'fasta')]
    lst_clip_record = list()
    for record in lst_seq_record:
        try:
            clip_rec = ClipRecord(record)
            lst_clip_record.append(clip_rec)
        except HeaderError:
            print('Sequence with invalid header')
            continue
    
    return lst_clip_record

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

def write_germ_tab_file(fname_tab,lst_seq_record):
    '''
    This function will take a list of Bio.SeqRecord instances meant to
    represent germline sequences, extract the needed sequence and metadata
    from each, and write them to an appropriately formatted tab file.
    '''
    with open(fname_tab,'w') as f:
        tab_writer = csv.DictWriter(f, tpl_germ_cols, extrasaction='ignore',
            delimiter='\t')
        # write header
        tab_writer.writeheader()
        for seq_record in lst_seq_record:
            tab_writer.writerow(seq_record.dict_germ_attr)
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

def mask_d_region_length(str_seq, v_length, j_length, mask_char='n'):
    '''
    This function masks the D region of a given sequence, V length, and J
    length.
    The optional argument 'mask_char' defaults to 'n'.
    '''
    length = len(str_seq)
    d_length = length - v_length - j_length
    lst_seq = list(str_seq)
    lst_seq[v_length:v_length + d_length] = d_length * mask_char
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
    return ''.join(lst_seq)

def parse_header_delim(str_line, field_sep='|', colon=':'):
    '''
    This function parses the header line of a FASTA file, scraping up its
    key-value pairs and returning them as a dictionary.
    It accepts 1 positional, and 2 optional args:
    str_line is the string form of the header line.
    field_sep is the separator in between key-value pairs (defaults to '|').
    colon is the separator between keys and values
    (defaults to ':', obviously).
    
    It expects header lines to be in the form of '[|]key:value[|key:value]'.
    '''
    # splits description by 'field_sep' into key-value pair strings
    # only noticing key-value strings longer than 0 characters
    lst_str_keyval = [keyval for keyval in
        str_line.split(field_sep) if len(keyval)>0]
    # splits each key value pair string by 'colon', and creates dictionary
    # NOTE: only split by the *first* 'colon' so that additional 'colon' can
    # appear in the value
    try:
        dict_header_field = dict([keyval.split(colon,1) for keyval in
            lst_str_keyval if len(keyval)>0])
    except:
        raise HeaderError(lst_str_keyval)
    return dict_header_field

def fasta2tab_file(fname, mask=None, germ=True):
    '''
    This function will process a single clip fasta file into a single tab file
    appropriately named.
    '''
    #generate an initial list of SeqRecord
    lst_seq_record = read_fasta_file(fname)
    #prune out entries that list germlines
    lst_germline, lst_clip_seq = prune_germline_records(lst_seq_record)
    #prep germline sequences
    dict_germline = dict()
    for germline in lst_germline:
        germline.prep_germ()
        dict_germline[germline.dict_germ_attr['clone']
            ] = germline.dict_germ_attr['seq']
    #parse supplementary attributes
    #mask d region (maybe)
    if mask=='same_null':
        for clip_seq in lst_clip_seq:
            clip_seq.d_no_mask()
    elif mask=='same_mask':
        for clip_seq in lst_clip_seq:
            clip_seq.d_mask()
    elif mask=='germ_null':
        for clip_seq in lst_clip_seq:
            clip_seq.germline_d_no_mask(dict_germline)
    else:
        for clip_seq in lst_clip_seq:
            clip_seq.d_no_mask()
    #append uuid portion to SEQUENCE_ID
    append_uuid(lst_clip_seq)
    #find an appropriate name for the tab file
    tabfname = tabname(fname)
    #write the final list of SeqRecord to the tab file
    write_tab_file(tabfname, lst_clip_seq)
    
    if germ:
        ## Only if preparing a germline tab
        #find an appropriate name for the germline tab file
        germtabfname = germtabname(fname)
        write_germ_tab_file(germtabfname, lst_germline)
    
    return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert FASTA files to TAB files'
        )
    parser.add_argument('files', metavar='infiles', nargs='+')
    parser.add_argument('-m', '--mask',
        choices=['same_null','same_mask','germ_null','germ_mask'])
    parser.add_argument('-g', '--germ', action='store_true')
    argspace = parser.parse_args()
    lst_fasta_files = argspace.files
    mask_mode = argspace.mask
    for fname in lst_fasta_files:
        fasta2tab_file(fname, mask=mask_mode)

