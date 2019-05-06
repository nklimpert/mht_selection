'can be run from within a paml directory'
from ete3 import EvolTree
import os

tree_file = "testTree.tre"
alignment_file = "testAlignment.fasta"
model = "./model/out"
model_name = "bsD.bl_0.2w"
# model_name = os.path.basename(os.getcwd())


testTree = EvolTree(tree_file)
testTree.link_to_alignment(alignment_file)
testTree.link_to_evol_model(model, model_name)


testTree.show()