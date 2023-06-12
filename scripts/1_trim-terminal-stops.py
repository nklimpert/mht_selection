'''Usage: trim-terminal-stops.py <alignment>'''

import os  # Manipulating filenames
from Bio import SeqIO


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
    alignment_file = 'atpE_petrosaviid_1_2018-09-29.fasta'
    
alignment_name = os.path.splitext(alignment_file)[0]  # Retrieve filename
alignment_outname = alignment_name + '.fasta' # + '_no-terminal-stops.fasta'
alignment_format = os.path.splitext(alignment_file)[1]  # Retrieve filetype
alignment_format = alignment_format[1:]  # Remove '.' character from filetype

terminal_stops_trimmed_alignment = []
for record in SeqIO.parse(alignment_file, format=alignment_format):
    record.seq = record.seq[:-3]
    terminal_stops_trimmed_alignment.append(record)


SeqIO.write(terminal_stops_trimmed_alignment, handle=alignment_outname, \
            format=alignment_format)