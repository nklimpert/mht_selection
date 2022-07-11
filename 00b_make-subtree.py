'''Usage: make-subtree.py <tree> <taxa_to_keep> <clade>'''

from ete3 import Tree  # Manipulating trees
import os  # Manipulating filenames

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
    tree_file = cmdln_args.get('<tree>')
    taxa_to_keep_file = cmdln_args.get('<taxa_to_keep>')
    clade = cmdln_args.get('<clade>')
# Run interactively in an iPython console
if in_ipython() is True:
    tree_file = 'petrosaviids_full.tre'
    taxa_to_keep_file = 'petrosaviids_1_taxa_corallorhizaStriat.txt'
    clade = 'corallorhizaStriat'

tree_name = os.path.splitext(tree_file)[0]
out_tree_filename = tree_name + '_' + clade + '.tre'


full_tree = Tree(tree_file)

taxa_to_keep = []
with open(taxa_to_keep_file, 'r') as taxa_list:
    for taxon in taxa_list:
        taxon = taxon.rstrip('\n')    
        taxa_to_keep.append(taxon)

full_tree.prune(taxa_to_keep, preserve_branch_length=True)

full_tree.unroot()

full_tree.write(outfile=out_tree_filename)
