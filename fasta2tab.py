#!/usr/bin/env python
import random
import string
import csv
import warnings
import os

from Bio import SeqIO

# the character set for sequence UUIDs
charset = string.letters + string.digits
# the number of characters to use for sequence UUIDs
n_char = 8
# these are the columns to output in the tab file
tpl_cols_0 = (
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
    'FASTA',
    )
tpl_cols_1 = (
    'order',
    'seqID',
    'functional',
    'in-frame',
    'stop',
    'mutation_invariate',
    'v_match',
    'v_length',
    'j_match',
    'j_length',
    'v_call',
    'j_call',
    'v_gap_length',
    'j_gap_length',
    'juncton_gap_length',
    'junction_nt',
    'junction_aa',
    'gap_method',
    'subject',
    'subset',
    'tissue',
    'disease',
    'date',
    'lab',
    'experimenter',
    'copy_number_close',
    'collapse_to_close',
    'copy_number_iden',
    'collapse_to_iden',
    'sequence',
    'germline',
    'cloneID',
        )
class HeaderError(ValueError):
    pass
class InvalidHeaderStyle(ValueError):
    pass
class ClipRecord(object):
    mask_char = 'n'
    field_sep = '|'
    colon = ':'
    header_rev = 1

    @classmethod
    def set_header_rev(cls, header_rev):
        if header_rev not in [0, 1]:
            raise InvalidHeaderStyle(
                    'Invalid header_rev: {:d}'.format(header_rev))
        else:
            setattr(cls, 'header_rev', header_rev)
        return None

    @classmethod
    def append_uuid(cls, lst_seq_record, **kwarg):
        '''
        Accepts a list of Bio.SeqRecord, compute a unique
        identifier for each, and append it to the 'SEQUENCE_ID' string of each
        SeqRecord.

        Parameters
        ----------
        length : int
            The length of UUID to append, default `n_char`.
        column : str
            The label of the column to which to append the UUID.
        '''
        n_id = len(lst_seq_record)
        lst_uuid = get_n_uuid(n_id, **kwarg)
        if cls.header_rev == 0:
            column = kwarg.pop('column', 'SEQUENCE_ID')
        elif cls.header_rev == 1:
            column = kwarg.pop('column', 'seqID')
        else:
            column = kwarg.pop('column', 'SEQUENCE_ID')
        for iI, record in enumerate(lst_seq_record):
            try:
                setattr(record, column,
                        getattr(record, column, '') + lst_uuid[iI])
            except:
                print(record.description)
                print(dir(record))
                raise

    def __init__(self, seq_record, **kwarg):
        header_rev = kwarg.pop('header_rev', False)
        if header_rev:
            self.header_rev = header_rev
        self.seq_record = seq_record
        self.description = self.seq_record.description
        self.seq = self.seq_record.seq
        self.SEQUENCE = str(self.seq_record.seq).lower()
        # set fields from FASTA header row
        self._set_from_header()
        if self.header_rev == 0:
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
        elif self.header_rev == 1:
            if hasattr(self, 'cloneID'):
                if self.cloneID[:9] == 'Germline:':
                    self.clone = self.cloneID[9:]
                else:
                    self.clone = self.cloneID
        else:
            raise InvalidHeaderStyle()
        return None

    def _set_from_header(self):
        '''
        This method scrapes input fields from the FASTA sequence description
        in the format of '[|]key:value[|key:value]', and adds them all as
        instance attributes.
        This method is run by __init__, is private, because it should never
        need to be run from outside the module.
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
        # if not has_attr(self,...)
        if self.header_rev == 0:
            v_length = int(self.V_LENGTH)
            j_length = int(self.J_LENGTH)
            self.GERMLINE_GAP_D_MASK = mask_d_region_length(
                self.SEQUENCE, v_length, j_length, mask_char=self.mask_char)
        elif self.header_rev == 1:
            v_length = int(self.v_length)
            j_length = int(self.j_length)
            self.germline_gap_d_mask = mask_d_region_length(
                self.sequence, v_length, j_length, mask_char=self.mask_char)
        else:
            raise InvalidHeaderStyle()
        pass

    def d_no_mask(self):
        '''
        This method takes the sequence data and copies it into
        GERMLINE_GAP_D_MASK without actually masking the D region.
        '''
        if self.header_rev == 0:
            self.GERMLINE_GAP_D_MASK = self.SEQUENCE
        elif self.header_rev == 1:
            self.germline_gap_d_mask = self.sequence
        else:
            raise InvalidHeaderStyle()
        return None
    def germline_d_no_mask(self, dict_germline):
        '''
        This method populates the GERMLINE_GAP_D_MASK field with the
        appropriate sequence from the germline entries.
        '''
        if self.header_rev == 0:
            try:
                self.GERMLINE_GAP_D_MASK = dict_germline[self.CLONE]
            except KeyError:
                print('Unable to find germline matching clone: {clone}'.format(
                    clone=self.CLONE))
                raise
        elif self.header_rev == 1:
            try:
                self.germline_gap_d_mask = dict_germline[self.clone]
            except KeyError:
                print('Unable to find germline matching clone: {clone}'.format(
                    clone=self.cloneID))
                raise
        else:
            raise InvalidHeaderStyle()
        pass

    def prep_germ(self):
        '''
        This method prepares a ClipRecord instance for entry into a germline
        TAB file.
        '''
        if self.header_rev == 0:
            if hasattr(self,'Germline'):
                self.SEQUENCE_ID = getattr(self,'Germline')
            elif hasattr(self,'>Germline'):
                self.SEQUENCE_ID = getattr(self,'>Germline')
            elif hasattr(self,'>GERMLINE'):
                self.SEQUENCE_ID = getattr(self,'>GERMLINE')
            else:
                print(self.description)
                warnings.warn(
                    '''The above header had no suitable 'Germline' field.''')
                # raise HeaderError('''ClipRecord object has no suitable 'Germline' field.''')
        elif self.header_rev == 1:
            if hasattr(self,'Germline'):
                self.seqID = getattr(self,'Germline')
            elif hasattr(self,'>Germline'):
                self.seqID = getattr(self,'>Germline')
            elif hasattr(self,'>GERMLINE'):
                self.seqID = getattr(self,'>GERMLINE')
            else:
                print(self.description)
                warnings.warn(
                    '''The above header had no suitable 'Germline' field.''')
                # raise HeaderError('''ClipRecord object has no suitable 'Germline' field.''')
        else:
            raise InvalidHeaderStyle()
        pass

    @classmethod
    def from_tabfile(cls, tabfname):
        '''
        Creates a list of ClipRecord instances from a pre-existing TAB file.
        '''
        with open(tabfname,'rU') as tabf:
            reader = csv.DictReader(tabf, delimiter='\t')
            lst_cliprec = list()
            for entrydict in reader:
                cliprec = cls.__new__(cls)
                for key, val  in entrydict.items():
                    try:
                        setattr(cliprec, key, val)
                    except TypeError:
                        print('Error setting attribute-------------------')
                        print(cliprec)
                        print(key)
                        print(val)
                        print('-------------------------------------------')
                        # raise
                        pass
                lst_cliprec.append(cliprec)
        return lst_cliprec
    pass

def prune_germline_records(lst_seq_record):
    lst_germline = [record for record in lst_seq_record if
        record.description[0] =='>']
    lst_non_germline = [record for record in lst_seq_record if
        record.description[0] !='>']
    return lst_germline, lst_non_germline


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
    [
    SEQUENCE_ID,
    SEQUENCE_GAP,
    GERMLINE_GAP_DMASK,
    CLONE_ID,
    ]
    or hopefully it will use the new-style headers.
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
    for record in lst_clip_record:
        setattr(record, 'fasta', os.path.basename(fname_fasta))

    return lst_clip_record

def _get_fields(lst_seq_record, **kwarg):
    """
    Get the fields that should always be included, plus the uppercase ones
    from a list of `Bio.SeqRecord` instances.

    Parameters
    ----------
    lst_seq_record : list
        A list of `Bio.SeqRecord` instances.
    header_rev : int
        Which revision of the header spec to use. 0 for old-style header
        spec, 1 for new-style header spec. (default: 1)

    Returns
    -------
    tpl_fields : tuple
        A tuple of the fields.
    """
    # Use header_rev to determine which set of headers to use.
    header_rev = kwarg.pop('header_rev', 1)
    if header_rev == 0:
        tpl_fields = list(tpl_cols_0)
    elif header_rev == 1:
        tpl_fields = list(tpl_cols_1)
    else:
        raise ValueError('Invlalid header_rev: {:d}'.format(header_rev))
    for seq_record in lst_seq_record:
        tpl_fields.extend([field for field in vars(seq_record)
            if (field not in tpl_fields)])
    return tpl_fields

def add_uuid_tabfile(tabinname, **kwarg):
    '''
    Append a UUID to each entry in a tabfile.
    '''
    ## If no taboutname specified, use the same file.
    ## This behavior may need to be changed later.
    taboutname = kwarg.pop('taboutname', tabinname)
    lst_cliprec = ClipRecord.from_tabfile(tabinname)
    ClipRecord.append_uuid(lst_cliprec)
    write_tab_file(taboutname, lst_cliprec)
    return None

def write_tab_file(fname_tab,lst_seq_record):
    '''
    Write a list of `Bio.SeqRecord` instances to a tabfile.

    Parameters
    ----------
    fname_tab : str
        A filename to which to write the tabfile.
    lst_seq_record : list
        A list of `Bio.SeqRecord` instances.

    Notes
    -----
    The fields that will be written to the tabfile are determined by the
    fields listed in `tpl_cols`, as well as those present in the
    `Bio.SeqRecord` instances that have all uppercase.
    '''
    tpl_fields = _get_fields(lst_seq_record)
    with open(fname_tab,'w') as f:
        tab_writer = csv.DictWriter(f, tpl_fields, extrasaction='ignore',
            delimiter='\t')
        # write header
        tab_writer.writeheader()
        for seq_record in lst_seq_record:
            tab_writer.writerow(vars(seq_record))
            pass

def write_germ_tab_file(fname_tab,lst_seq_record):
    '''
    This function will take a list of Bio.SeqRecord instances meant to
    represent germline sequences, extract the needed sequence and metadata
    from each, and write them to an appropriately formatted tab file.
    '''
    with open(fname_tab,'w') as f:
        tab_writer = csv.DictWriter(f, tpl_cols, extrasaction='ignore',
            delimiter='\t')
        # write header
        tab_writer.writeheader()
        for seq_record in lst_seq_record:
            tab_writer.writerow(vars(seq_record))
            pass

def get_n_uuid(n, **kwarg):
    '''
    This function will generate as many UUID as requested, confirming that
    they are unique.

    Parameters
    ----------
    length : int
        The length of UUID to append, default `n_char`.
    '''
    def gen_uuid():
        return get_uuid(**kwarg)
    # generates list of UUID prior to checking for duplicates
    lst_uuid = [gen_uuid() for iI in xrange(n)]
    # WARNING: while loop, remove in future versions
    # check to see that all UUID in lst_uuid are unique
    while True:
        lst_uuid = set(lst_uuid)
        n_unique = len(lst_uuid)
        if n_unique == n:
            break
        lst_uuid.update([gen_uuid() for iI in xrange(n-n_unique)])
    return list(lst_uuid)

def get_uuid(**kwarg):
    '''
    This function returns a single UUID using module variables charset,n_char.

    Parameters
    ----------
    length : int
        The length of UUID to append, default `n_char`.
    '''
    length = kwarg.pop('length', n_char)
    return ''.join(random.choice(charset) for iI in xrange(length))

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
            # raise ValueError('''The sequence has changed in length.
            #     You must mask a different way.''')
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
        lst_tpl_keyval = [keyval.split(colon,1) for keyval in
            lst_str_keyval if len(keyval)>0]
        if len(lst_tpl_keyval) == 1 and len(lst_tpl_keyval[0]) == 1:
            # this means we are dealing with a thin header
            dict_header_field = {'SEQUENCE_ID':lst_tpl_keyval[0][0]}
        else:
            dict_header_field = dict(lst_tpl_keyval)
    except:
        raise HeaderError(lst_str_keyval)
    return dict_header_field

def fasta2tab_comb(lst_fname, tabfname='tabfile.tab', mask=None,
        germ=True, **kwarg):
    """
    Combine multiple clip fasta files `lst_fname` into a single tab file.
    """
    lst_seq_record = list()
    for fname in lst_fname:
        lst_seq_record.extend(read_fasta_file(fname))
    #process lst_seq_record
    seqrecords2tab(lst_seq_record, mask=mask, germ=germ,
            tabfname=tabfname, **kwarg)
    return None

def fasta2tab_file(fname, mask=None, germ=True, **kwarg):
    '''
    This function will process a single clip fasta file into a single tab file
    appropriately named.
    '''
    #generate an initial list of SeqRecord
    lst_seq_record = read_fasta_file(fname)
    #find an appropriate name for the tab file
    tabfname = tabname(fname)
    #process lst_seq_record
    seqrecords2tab(lst_seq_record, mask=mask, germ=germ,
            tabfname=tabfname, **kwarg)
    return None

def seqrecords2tab(lst_seq_record, **kwarg):
    """
    Process a list of `ClipRecord` instances into a tabfile.

    Parameters
    ----------
    lst_seq_record : list
        A list of `ClipRecord` instances, including germlines and
        non-germlines which will be processed.
    mask : str
        A string specifying what masking protocol to use. (default: None)
    germ : bool
        Whether or not to prepare a germline tab file. (default: True)
    integrate : bool
        Whether or not to integrate the germline sequences into the
        master tab file. (default: True)

    Returns
    -------
    None
    """
    #parse kwarg
    mask = kwarg.pop('mask', None)
    germ = kwarg.pop('germ', True)
    tabfname = kwarg.pop('tabfname', True)
    integrate = kwarg.pop('integrate', True)
    #prune out entries that list germlines
    lst_germline, lst_clip_seq = prune_germline_records(lst_seq_record)
    #prep germline sequences
    dict_germline = dict()
    for germline in lst_germline:
        germline.prep_germ()
        dict_germline[getattr(germline,'SEQUENCE_ID','anonymous')
            ] = getattr(germline, 'SEQUENCE')
    #parse supplementary attributes
    mask_by_option(lst_clip_seq, dict_germline=dict_germline, mask=mask)
    #append uuid portion to SEQUENCE_ID
    ClipRecord.append_uuid(lst_clip_seq)
    if integrate:
        # If told to integrate germline sequences, prepend those to the
        # list of clip seq.
        lst_clip_seq = lst_germline + lst_clip_seq
    #write the final list of SeqRecord to the tab file
    write_tab_file(tabfname, lst_clip_seq)

    if germ:
        ## Only if preparing a germline tab
        #find an appropriate name for the germline tab file
        germtabfname = germtabname(tabfname)
        write_germ_tab_file(germtabfname, lst_germline)

    return None

def mask_by_option(lst_clip_seq, dict_germline=None, mask=None):
    """
    Mask sequences to fill 'GERMLINE_GAP_D_MASK' field.

    Parameters
    ----------

    lst_clip_seq : list
        List of `ClipRecord` instances to mask by option

    Returns
    -------
    lst_clip_seq : list
        List of `ClipRecord` instances that have been masked.
    """
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
        # same as if mask=='same_null'
        for clip_seq in lst_clip_seq:
            clip_seq.d_no_mask()
    return lst_clip_seq

def _main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert FASTA files to TAB files'
        )
    parser.add_argument('files', metavar='infiles', nargs='+',
            help="""
            FASTA files to convert into TAB files.
            They can be either treated individually or as a group,
            depending on '-c'.
            """,
            )
    parser.add_argument('-m', '--mask',
            choices=['same_null','same_mask','germ_null','germ_mask'],
            default='same_null',
            help="""
            DEPRECATED:
            This option is to modify how the column
            'GERMLINE_GAP_D_MASK' is populated.
            There is no reason to use this option unless you are
            generating a TAB file to work with the legacy auto phylip
            program written in R by Jason Van Der Heiden.
            (default: 'same_null')
            """,
            )
    parser.add_argument('-g', '--germ',
            dest='germ',
            action='store_true',
            help="""
            Whether or not the generate a separate TAB file that
            includes germline sequences found in the input FASTA files.
            This option was required for an older workflow, where
            germline sequences would be manually integrated into output
            TAB files. (default: False)
            """,
            )
    parser.add_argument('-c', '--combine',
            dest='tabfname',
            nargs='?', const='tabfile.tab', default=None,
            help="""
            If specified, it treats the input FASTA files as
            meaningfully related, and joints them into a single TAB
            file.
            If a filename is provided as well, that name is used,
            otherwise a default name is used. (default: 'tabfile.tab')
            """,
            )
    parser.add_argument('-i', '--integrate',
            dest='integrate',
            action='store_true',
            help="""
            Whether or not to integrate germline sequences from the
            input FASTA files into the output TAB file[s].
            Germline sequences are output all prior to regular
            sequences, and they do not receive UUIDs, so they can have
            collisions depending on their names.
            The way this option works may be changed in future versions.
            (default: False)
            """,
            )
    parser.add_argument('-r', '--header',
            dest='header',
            default=1,
            type=int,
            help="""
            If specified, one can change the version of headers to use.
            Old-style headers are 0, new-style headers are 1 (default).
            """,
            )
    argspace = parser.parse_args()
    mask_mode = argspace.mask
    ClipRecord.set_header_rev(argspace.header)
    if argspace.tabfname == None:
        # No output tabfile name was specified, therefore treat each input
        # fasta file individually.
        for fname in argspace.files:
            fasta2tab_file(
                fname,
                mask=argspace.mask,
                germ=argspace.germ,
                integrate=argspace.integrate,
                )
    else:
        fasta2tab_comb(
            argspace.files,
            tabfname=argspace.tabfname,
            mask=argspace.mask,
            germ=argspace.germ,
            integrate=argspace.integrate,
            )

def _tabmod():
    import argparse
    parser = argparse.ArgumentParser(
        description='Modify TAB files in various ways'
        )
    parser.add_argument('lst_file', metavar='FILES',
            nargs='+',
            help="""The TAB file(s) to which to append UUIDs.
            """
            )
    parser.add_argument('-r', '--header',
            dest='header',
            default=1,
            type=int,
            help="""
            If specified, one can change the version of headers to use.
            Old-style headers are 0, new-style headers are 1 (default).
            """,
            )
    parser.add_argument('-o', '--output',
            dest='output',
            help="""
            If specified, place the modified TAB file at this path
            (instead of in place).
            This option should only be used for 1 at a time TAB files.
            """,
            )
    subparsers = parser.add_subparsers(help="""
            Choose from available subcommands
            """,
            dest='action',
            )
    parser_uuid = subparsers.add_parser('u', help="""
            Add UUIDs to SEQUENCE_ID field of TAB file(s).
            """,
            )
    parser_uuid.add_argument('-l', '--length',
            default=n_char,
            dest='length',
            help="""The number of characters for each UUID to have
            (default {num}).
            The character set from which the UUID is chosen is: '{chars}'.
            """.format(num=n_char,chars=charset),
            )
    # parser_uuid.set_defaults(func=add_uuid_tabfile)
    argspace = parser.parse_args()
    if argspace.action == 'u':
        ClipRecord.set_header_rev(argspace.header)
        length = getattr(argspace, 'length', n_char)
        if hasattr(argspace, 'output') and len(argspace.lst_file) > 1:
            print('''Output file should not be specified with more than
            one input file.''')
            raise
        if len(argspace.lst_file) == 1 and argspace.output:
            tabfile = argspace.lst_file[0]
            add_uuid_tabfile(tabfile,
                    taboutname=getattr(argspace, 'output', tabfile),
                    length=length,
                    )
        else:
            for tabfile in argspace.lst_file:
                add_uuid_tabfile(tabfile,
                        length=length,
                        )
    return None

if __name__ == '__main__':
    _main()
