from Bio import AlignIO
import os

missing = {}

for alignment in [s for s in os.listdir(os.getcwd()) if s.endswith(".fasta")]:
    gene = alignment.split("_")[0]

    # loop through all of the records
    for record in AlignIO.read(alignment, "fasta"):
        if(str(record.seq) == '-'*len(record.seq)): # if it's just gaps
            if(record.name not in missing):
                missing[record.name] = [gene]
            else:
                missing[record.name].append(gene)





# with open("missingGreens.txt", 'w') as outFile:
#     for key in missing.keys():
#         outFile.write(key+"\n")
#         for value in missing[key]:
#             outFile.write(value+ "\n")
