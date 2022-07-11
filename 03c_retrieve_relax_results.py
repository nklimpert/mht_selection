'''Usage: retrieve_results.py <json_file>'''

# Modules
import json


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
    json_file = cmdln_args.get('<json_file>')
# Run interactively in an iPython console
if in_ipython() is True:
    json_file = 'accD_gymnoApter_2018-10-02.fasta.RELAX.json'
    

alignment_name = json_file.split('.')[0]
gene_name = alignment_name.split('_')[0]
clade_name = alignment_name.split('_')[1]
date = alignment_name.split('_')[2]

with open(json_file, "r") as JSON:
    results = json.load(JSON)

LRT_results = results.get('relaxation-test')
p_value = LRT_results.get('p')

fits_results = results.get('fits')
alternative_results = fits_results.get('Alternative')
K_value = alternative_results.get('K')

final_results = clade_name + ',' + gene_name + ',' + str(p_value) + ',' + \
                str(K_value) + '\n'

with open('results.csv', 'w') as results:
    results.write(final_results)
