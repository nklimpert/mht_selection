'''Usage: make_relax_tree.py <tree> <alignment> <test_taxa>'''

from Bio import SeqIO
from ete3 import EvolTree  # Manipulating trees
import os  # Manipulating filenames
from itertools import combinations

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
    alignment_file = cmdln_args.get('<alignment>')
    test_taxa_file = cmdln_args.get('<test_taxa>')
# Run interactively in an iPython console
if in_ipython() is True:
    tree_file = 'petrosaviids_aphyllorchis.tre'
    alignment_file = 'accD_aphyllorchis_2018-10-02.fasta'
    test_taxa_file = 'test_taxa_aphyllorchis.txt'


alignment_name = os.path.basename(alignment_file)  # Retrieve filename
gene_name = alignment_name.split('_')[0]  # Retrieve gene name
alignment_format = os.path.splitext(alignment_file)[1]  # Retrieve filetype
alignment_format = alignment_format[1:]  # Remove '.' character from filetype


# Check if the alignment contains any empty sequences
empty_seq_count = 0
for record in SeqIO.parse(alignment_file, format=alignment_format):
    gapSeq = '-'*len(record.seq)
    if (str(record.seq).upper().replace("N", "-")) == gapSeq:  # if it's just gaps
        empty_seq_count += 1

# If there were empty sequences found in the alignment, make a new alignment
# with empty sequences removed and record the names of the taxa in the new 
# alignment for pruning the tree
if empty_seq_count >= 1:
    trimmed_alignment = []
    taxa_in_alignment = []
    for record in SeqIO.parse(alignment_file, format=alignment_format):
        gapSeq = '-'*len(record.seq)
        if (str(record.seq).upper().replace("N", "-")) == gapSeq:  # if it's just gaps
            pass
        else:
            trimmed_alignment.append(record)
            taxa_in_alignment.append(record.id)

# Only write a new alignment if there is a new alignment
if empty_seq_count >= 1:
    if len(trimmed_alignment) >= 1:
        SeqIO.write(trimmed_alignment, handle=alignment_file, \
                    format=alignment_format)


tree = EvolTree(tree_file)
out_tree_name = os.path.basename(tree_file)
out_tree_name = os.path.splitext(out_tree_name)[0]
out_tree_name = out_tree_name + '_' + gene_name + '.tre'

# If there is a new alignment, prune the tree down to the taxa that remain in
# the new alignment
if empty_seq_count >= 1:
    if len(taxa_in_alignment) >= 1:
        tree.prune(taxa_in_alignment, preserve_branch_length=True)

test_taxa = []
with open(test_taxa_file, 'r') as test_taxa_list:
    for taxon in test_taxa_list:
        taxon = taxon.rstrip()
        test_taxa.append(taxon)


nodes_to_mark = set() # set since we want it to be all unique ids

# Mark the test taxa
for taxon in test_taxa:
    taxon_node = tree&taxon # ete3 notation for finding a node within a tree
    taxon_id = taxon_node.node_id # mark_tree only takes node_ids, not labels
    nodes_to_mark.add(taxon_id)

# Find internal nodes below the test taxa and mark them
for i in range(len(test_taxa), 1,-1):
    taxa_groups = [x for x in combinations(test_taxa, i)]
    for group in taxa_groups:
        common_node = tree.get_common_ancestor(*group)
        taxon_id = common_node.node_id
        nodes_to_mark.add(taxon_id)

#TODO change the names of the nodes
for mark_id in nodes_to_mark:
    test_node = tree.search_nodes(node_id=mark_id)[0]
    test_node.name += '{test}'

tree.write(outfile=out_tree_name, format=1)
