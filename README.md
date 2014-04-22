# fasta2tab

**Jason Gilliland**

This is a utility program to help convert a CLIP FASTA file into a TAB
file for use with some other sequence analysis programs.

Usage
-----
```
usage: fasta2tab.py [-h] [-m {same_null,same_mask,germ_null,germ_mask}] [-g]
                    [-c [TABFNAME]]
                    infiles [infiles ...]

Convert FASTA files to TAB files

positional arguments:
  infiles

optional arguments:
  -h, --help            show this help message and exit
  -m {same_null,same_mask,germ_null,germ_mask}, --mask {same_null,same_mask,germ_null,germ_mask}
  -g, --germ
  -c [TABFNAME], --combine [TABFNAME]
```

Each file that you provide to fasta2tab.py should be a fasta (\*.fasta) file
with a number of esoteric key-value pairs in the description.
For each file that you provide, fasta2tab.py will generate a tab delimited
values file (\*.tab) which will place some of those keys in appropriate
columns, as well as compute some new keys from some of the fields provided.

It depends on:
- BioPython (a dependency which may be removed at a later date)
