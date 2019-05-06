'''Usage: check_frame.py <alignment>'''

import os  # Manipulating filenames
from Bio import SeqIO
from Bio.Data import CodonTable

# Check if running interactively in an iPython console, or in a script from the
# command line
def in_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
# Run in a script from the command line
if in_ipython() is False:
    from docopt import docopt  # Command line argument handler
    cmdln_args = docopt(__doc__)
    alignment_file = cmdln_args.get('<alignment>')
# Run interactively in an iPython console
if in_ipython() is True:
    alignment_file = 'rpl20_petrosavia_2018-10-01.fasta'

alignment_format = os.path.splitext(alignment_file)[1]  # Retrieve filetype
alignment_format = alignment_format[1:]  # Remove '.' character from filetype

record = next(SeqIO.parse(alignment_file, format=alignment_format))
seq_len = len(record.seq)
if seq_len % 3 != 0:
    print alignment_file + ' is out of frame.'
else:
    for record in SeqIO.parse(alignment_file, format=alignment_format):
        try:
            aa_record = record.translate(gap='-')
            if '*' in aa_record.seq:
                print alignment_file + ' contains stop codons.'
        except CodonTable.TranslationError:
            pass
            print alignment_file + ' has a Translation Error.'


#record_nogaps = record.seq.ungap('-')
#if len(record_nogaps) % 3 != 0:
#   pass
#else:
#   aa_record_nogaps = record_nogaps.translate()
#if '*' in aa_record_nogaps:
#print alignment_file + ' contains stop codons.'
#else:
#if len(record.seq) % 3 !=0:
#pass
#else: