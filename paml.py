'''Usage: paml.py <alignment> <tree> <test> <test_taxa>'''

from Bio import SeqIO
import os
from ete3 import EvolTree
from itertools import combinations

tests = {'branch': ['M0', 'b_free'],
         'bsA': ['bsA1', 'bsA'],
         'cmD': ['M3', 'bsD'],
         'cmC': ['XX', 'bsC']
         } # NOTE XX refers to M2a_rel - this is just how to do a user-defined model

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
    tree_file = cmdln_args.get('<tree>')
    test = tests[cmdln_args.get('<test>')]
    test_taxa_file = cmdln_args.get('<test_taxa>')
# Run interactively in an iPython console
if in_ipython() is True:
    alignment_file = 'infA_aphyllorchis_2018-10-29.fasta'
    tree_file = 'petrosaviids_aphyllorchis.tre'
    test = tests['branch']
    test_taxa_file = 'test_taxa_aphyllorchis.txt'

alignment_name = os.path.basename(alignment_file)  # Retrieve filename
gene_name = alignment_name.split('_')[0]  # Retrieve gene name
clade_name = alignment_file.split('_')[1]
alignment_format = os.path.splitext(alignment_file)[1]  # Retrieve filetype
alignment_format = alignment_format[1:]  # Remove '.' character from filetype


# Check if the alignment contains any empty sequences
empty_seq_count = 0
for record in SeqIO.parse(alignment_file, format=alignment_format):
    gapSeq = '-'*len(record.seq)
    if (str(record.seq).upper().replace("N", "-")) == gapSeq:  # if it's just gaps
        empty_seq_count += 1

# If there were empty sequences found in the alignment, record the names of the
# taxa with sequences for pruning the tree
if empty_seq_count >= 1:
    taxa_in_alignment = []
    for record in SeqIO.parse(alignment_file, format=alignment_format):
        gapSeq = '-'*len(record.seq)
        if (str(record.seq).upper().replace("N", "-")) == gapSeq:  # if it's just gaps
            pass
        else:
            taxa_in_alignment.append(record.id)

tree = EvolTree(tree_file)
out_tree_name = os.path.basename(tree_file)
out_tree_name = os.path.splitext(out_tree_name)[0]
out_tree_name = out_tree_name + '_' + gene_name + '.tre'

# If there is a new alignment, prune the tree down to the taxa that remain in
# the new alignment and write a new tree because EvolTree is shit and can't 
# use the pruned tree saved in memory
if empty_seq_count >= 1:
    if len(taxa_in_alignment) >= 1:
        tree.prune(taxa_in_alignment, preserve_branch_length=True)
        tree.unroot()
        tree.write(outfile=out_tree_name, format=0)
        tree = EvolTree(out_tree_name)

tree.link_to_alignment(alignment_file)

tree.workdir = os.getcwd()

# Record list of all node_ids in the tree for later retrieving omega from a
# background branch in the b_free model
list_of_node_ids = []
for node in tree.traverse('postorder'):
    list_of_node_ids.append(node.node_id)

test_taxa = []
with open(test_taxa_file, 'r') as test_taxa_list:
    for taxon in test_taxa_list:
        taxon = taxon.rstrip()
        test_taxa.append(taxon)

marked_taxon_ids = []
# Mark test taxa
for taxon in test_taxa:
    taxon_node = tree&taxon # ete3 notation for finding a node within a tree
    marked_taxon_id = taxon_node.node_id # mark_tree only takes node_ids, not 
                                         # labels
    tree.mark_tree([marked_taxon_id])
    marked_taxon_ids.append(marked_taxon_id)

# Find internal nodes below the test taxa and mark them
for i in range(len(test_taxa), 1,-1):
    taxa_groups = [x for x in combinations(test_taxa, i)]
    for group in taxa_groups:
        common_node = tree.get_common_ancestor(*group)
        marked_taxon_id = common_node.node_id
        tree.mark_tree([marked_taxon_id])
        marked_taxon_ids.append(marked_taxon_id)

best_model = {'M0': None, 'b_free': None, 'bsA1': None, 'bsA': None, 
              'M3': None, 'bsD': None, 'XX': None, 'bsC': None}
best_lnL = {'M0': float('-inf'), 'b_free': float('-inf'), 
            'bsA1': float('-inf'), 'bsA': float('-inf'), 'M3': float('-inf'), 
            'bsD': float('-inf'), 'XX': float('-inf'), 
            'bsC': float('-inf')}

# Quicker version of running PAML for testing
# for model in test:
#     model_specifications = model + '.' + 'bl' + '_' + '0.7' + 'w'
#     print 'Testing model ' + model + ' on ' + alignment_name + \
#           ' using starting branch length option ' + \
#           'bl' + ' and initial omega: ' + \
#           '0.7' + 'w'
#     tree.run_model(model_specifications, fix_blength=1, omega=0.7)
#     current_model = tree.get_evol_model(model_specifications)
#     print 'The fitting of model ' + model + ' on ' + alignment_name + \
#            ' using starting branch length option ' + \
#            'bl' + ' and initial omega: ' + \
#            '0.7' + 'w, the likelihood was: ' + \
#            str(current_model.lnL)
#     if current_model.lnL > best_lnL[model]:
#         best_lnL[model] = current_model.lnL
#         best_model[model] = current_model

# Loop for running each version of each model in one test. Best models are
# stored in a dictionary that is later read to write results
for model in test:
    for starting_branch_length_option in [1, -1]: 
        if starting_branch_length_option == 1:
            branch_estimation = 'bl'
        elif starting_branch_length_option == -1:
            branch_estimation = 'random'
        for initial_omega in [0.2, 0.7, 1.2]:
            if model == 'bsA1':
                initial_omega = 1.0
            model_specifications = model + '.' + branch_estimation + '_' + \
                                   str(initial_omega) + 'w'
            print('Testing model ' + model + ' on ' + alignment_name + \
                  ' using starting branch length option ' + \
                  branch_estimation + ' and initial omega: ' + \
                  str(initial_omega) + 'w')
            if model == 'XX':
                tree.run_model(model_specifications, \
                            fix_blength=starting_branch_length_option, \
                            omega=initial_omega, NSsites=22, ncatG=3)
                            
                            
                # Here's the garbage I wrote to make sure that it parses the out files correctly
                tree.get_evol_model(model_specifications).properties['typ'] = 'branch-site'
                tree.get_evol_model(model_specifications)._load(model_specifications+'/out')
                
                
            else:
                tree.run_model(model_specifications, \
                            fix_blength=starting_branch_length_option, \
                            omega=initial_omega)
            current_model = tree.get_evol_model(model_specifications)
            print('The fitting of model ' + model + ' on ' + alignment_name + \
                  ' using starting branch length option ' + \
                  branch_estimation + ' and initial omega: ' + \
                  str(initial_omega) + 'w, the likelihood was: ' + \
                  str(current_model.lnL))
            if current_model.lnL > best_lnL[model]:
                best_lnL[model] = current_model.lnL
                best_model[model] = current_model
            if model == 'bsA1':
                break

# Retrieve results from the best model
for model in test:
    current_model = best_model[model]
    model_name = current_model.name
    lnL = current_model.lnL
    if model == 'M0':
        all_branch_stats = current_model.branches
        one_branch = all_branch_stats[1]
        omega = one_branch.get('w')
        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(omega) + '\n' 
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'b_free':
        all_branch_stats = current_model.branches
        fg_branch = all_branch_stats[marked_taxon_id]
        if fg_branch.get('mark') == ' #1':
            fg_omega = fg_branch.get('w')
        for node_id in list_of_node_ids:
            if node_id not in marked_taxon_ids:
                bg_branch_id = node_id
                bg_branch = current_model.branches[bg_branch_id]
                if bg_branch.get('mark') == ' #0':
                    bg_omega = bg_branch.get('w')
                    break
        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(bg_omega) + ',' + str(fg_omega) + '\n' 
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'M3':
        classes = current_model.classes
        proportions = classes.get('proportions')
        omegas = classes.get('w')
        proportion_0 = proportions[0]
        omega_0 = omegas[0]
        proportion_1 = proportions[1]
        omega_1 = omegas[1]
        proportion_2 = proportions[2]
        omega_2 = omegas[2]
        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(proportion_0) + ',' + str(omega_0) + \
                  ',' + str(proportion_1) + ',' + str(omega_1) + ',' + \
                  str(proportion_2) + ',' + str(omega_2) + '\n' 
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'bsD':
        classes = current_model.classes
        proportions = classes.get('proportions')
        background_omegas = classes.get('branch type 0')
        foreground_omegas = classes.get('branch type 1')
        proportion_0 = proportions[0]
        bg_omega_0 = background_omegas[0]
        fg_omega_0 = foreground_omegas[0]
        proportion_1 = proportions[1]
        bg_omega_1 = background_omegas[1]
        fg_omega_1 = foreground_omegas[1]
        proportion_2 = proportions[2]
        bg_omega_2 = background_omegas[2]
        fg_omega_2 = foreground_omegas[2]
        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(proportion_0) + ',' + str(bg_omega_0) +\
                  ',' + str(fg_omega_0) + ',' + str(proportion_1) + ',' + \
                  str(bg_omega_1) + ',' + str(fg_omega_1) + ',' + \
                  str(proportion_2) + ',' + str(bg_omega_2) + ',' + \
                  str(fg_omega_2) + '\n' 
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'XX': #If model is M2a_rel
        classes = current_model.classes
        proportions = classes.get('proportions')
        omegas = classes.get('w')
        proportion_0 = proportions[0]
        omega_0 = omegas[0]
        proportion_1 = proportions[1]
        omega_1 = omegas[1]
        proportion_2 = proportions[2]
        omega_2 = omegas[2]
        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(proportion_0) + ',' + str(omega_0) + \
                  ',' + str(proportion_1) + ',' + str(omega_1) + ',' + \
                  str(proportion_2) + ',' + str(omega_2) + '\n' 
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'bsC':
        classes = current_model.classes
        proportions = classes.get('proportions')
        background_omegas = classes.get('branch type 0')
        foreground_omegas = classes.get('branch type 1')
        proportion_0 = proportions[0]
        bg_omega_0 = background_omegas[0]
        fg_omega_0 = foreground_omegas[0]
        proportion_1 = proportions[1]
        bg_omega_1 = background_omegas[1]
        fg_omega_1 = foreground_omegas[1]
        proportion_2 = proportions[2]
        bg_omega_2 = background_omegas[2]
        fg_omega_2 = foreground_omegas[2]
        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(proportion_0) + ',' + str(bg_omega_0) +\
                  ',' + str(fg_omega_0) + ',' + str(proportion_1) + ',' + \
                  str(bg_omega_1) + ',' + str(fg_omega_1) + ',' + \
                  str(proportion_2) + ',' + str(bg_omega_2) + ',' + \
                  str(fg_omega_2) + '\n' 
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'bsA':
        classes = current_model.classes
        proportions = classes.get('proportions')
        background_omegas = classes.get('background w')
        foreground_omegas = classes.get('foreground w')
        proportion_0 = proportions[0]
        bg_omega_0 = background_omegas[0]
        fg_omega_0 = foreground_omegas[0]
        proportion_1 = proportions[1]
        bg_omega_1 = background_omegas[1]
        fg_omega_1 = foreground_omegas[1]
        proportion_2 = proportions[2]
        bg_omega_2a = background_omegas[2]
        fg_omega_2a = foreground_omegas[2]
        bg_omega_2b = background_omegas[3]
        fg_omega_2b = foreground_omegas[3]

        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(proportion_0) + ',' + str(bg_omega_0) +\
                  ',' + str(fg_omega_0) + ',' + str(proportion_1) + ',' + \
                  str(bg_omega_1) + ',' + str(fg_omega_1) + ',' + \
                  str(proportion_2) + ',' + str(bg_omega_2a) + ',' + \
                  str(fg_omega_2a) + str(bg_omega_2b) + ',' + \
                  str(fg_omega_2b) + '\n'
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
    if model == 'bsA1':
        classes = current_model.classes
        proportions = classes.get('proportions')
        background_omegas = classes.get('background w')
        foreground_omegas = classes.get('foreground w')
        proportion_0 = proportions[0]
        bg_omega_0 = background_omegas[0]
        fg_omega_0 = foreground_omegas[0]
        proportion_1 = proportions[1]
        bg_omega_1 = background_omegas[1]
        fg_omega_1 = foreground_omegas[1]
        proportion_2 = proportions[2]
        bg_omega_2a = background_omegas[2]
        fg_omega_2a = foreground_omegas[2]
        bg_omega_2b = background_omegas[3]
        fg_omega_2b = foreground_omegas[3]

        results = clade_name + ',' + gene_name + ',' + model_name + ',' + \
                  str(lnL) + ',' + str(proportion_0) + ',' + str(bg_omega_0) +\
                  ',' + str(fg_omega_0) + ',' + str(proportion_1) + ',' + \
                  str(bg_omega_1) + ',' + str(fg_omega_1) + ',' + \
                  str(proportion_2) + ',' + str(bg_omega_2a) + ',' + \
                  str(fg_omega_2a) + str(bg_omega_2b) + ',' + \
                  str(fg_omega_2b) + '\n'
        out_filename = clade_name + '_' + gene_name + '_' + model + '.csv'
        with open(out_filename, 'w') as out_results:
            out_results.write(results)
