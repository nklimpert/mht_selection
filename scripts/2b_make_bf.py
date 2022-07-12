'''Usage: make_bf.py <alignment>'''


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
    alignment_file = 'cemA_petrosaviid_1_2018-09-29.fasta'

gene_name = alignment_file.split('_')[0]  # Retrieve gene name
lineage_name = alignment_file.split('_')[1]

start_string = '''fileToExe = "/usr/lib/hyphy/TemplateBatchFiles/RELAX.bf";
inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]= "/home/wesley/ownCloud/research/mycoheterotroph-selection/relax/'''

middle_string = lineage_name + '/' + gene_name + '/' + alignment_file

end_string = '''";
inputRedirect["03"]="Y";
inputRedirect["04"]="test";
inputRedirect["05"]="Minimal";
ExecuteAFile( fileToExe, inputRedirect);'''

final_string = start_string + middle_string + end_string

with open('relax_run.bf', 'w') as outfile:
    outfile.write(final_string)
