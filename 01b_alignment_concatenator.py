from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
from sys import argv
from glob import glob
from Bio.Seq import UnknownSeq


def check_gene_group(geneName, geneDict='plastid'):
    '''takes a gene name and looks it up against a dictionary of gene names separated into groups.
       if no gene dictionary is given, it will use a preset dictionary.'''
    if geneDict == 'plastid':
        geneDict = {'accD': ['accD'], 'atp': ['atpA', 'atpB', 'atpE', 'atpF', 'atpH', 'atpI'],
                    'ccsA': ['ccsA'], 'cemA': ['cemA'], 'clpP': ['clpP'], 'infA': ['infA'], 'matK': ['matK'], 'misc': ['accD', 'matK', 'clpP'],
                    'ndh': ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK'],
                    'pet': ['petA', 'petB', 'petD', 'petG', 'petL', 'petN'], 'psa': ['psaA', 'psaB', 'psaC', 'psaI', 'psaJ', 'psbA'],
                    'psb': ['psbB', 'psbC', 'psbD', 'psbE', 'psbF', 'psbH', 'psbI', 'psbJ', 'psbK', 'psbL', 'psbM', 'psbN', 'psbT', 'lhbA'],
                    'rbcL': ['rbcL'], 'rpl': ['rpl14', 'rpl16', 'rpl2', 'rpl20', 'rpl22', 'rpl23', 'rpl32', 'rpl33', 'rpl36'],
                    'rpo': ['rpoA', 'rpoB', 'rpoC1', 'rpoC2'], 'rps': ['rps11', 'rps12', 'rps14', 'rps15', 'rps16', 'rps18', 'rps19', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8'],
                    'ycf3': ['ycf3'], 'ycf4': ['ycf4']
                    }
    else:
        geneDict = {'mito': ['atp1', 'atp4', 'atp6', 'atp8', 'atp9', 'ccmb', 'ccmC', 'ccmFc', 'ccmFn', 'cob', 'cox1', 'cox2', 'cox3', 'matR',
                             'mttB', 'nad1', 'nad2', 'nad3', 'nad4L', 'nad4', 'nad5', 'nad6', 'nad7', 'nad9', 'rpl16', 'rpl2', 'rpl5', 'rps12',
                              'rps19', 'rps1', 'rps3', 'rps4', 'rps7', 'sdh4', ]}
    geneGroups = []
    for group, genes in geneDict.items():
        if geneName in genes:
            geneGroups.append(group)
    return geneGroups


def get_gene_name(alignmentFile):
    alignmentName = os.path.basename(alignmentFile)  # Retrieve filename
    geneName = alignmentName.split('_')[0]  # Retrieve gene name
    return geneName


def check_missing(alignment):
    '''alignment is a MultipleSeqAlignment object.
        returns a list of taxa that don't have a sequence in the alignment'''
    missing_taxa = []
    for record in alignment:
        gapSeq = '-'*len(alignment[0])
        if (str(record.seq).upper().replace("N", "-")) == gapSeq:
            missing_taxa.append(record.id)
    return missing_taxa


def remove_missing(alignment, taxaToRemove):
    '''alignment is a MultipleSeqAlignment object.
        taxaToRemove is a list of the ids of taxa to remove.'''
    indices = [i for i in range(len(alignment)) if alignment[i].id in taxaToRemove]
    return MultipleSeqAlignment([alignment[i] for i in range(len(alignment)) if i not in indices])


def replace_missing(alignment):
    '''alignment is a MultipleseqAlignment object.
       taxaToRemove is a list of the ids of taxa to replace with N's'''
    for record in alignment:
        gapSeq = '-'*len(alignment[0])
        if (str(record.seq).upper().replace("N", "-")) == gapSeq:
            record.seq = UnknownSeq(len(record.seq), character='N')
    return alignment


def get_alignment_type(alignmentFile):
    '''returns the file type of the alignment so that AlignIO can read it correctly'''
    ext = os.path.splitext(alignmentFile)[1].replace('.', '')
    return ext


def concat(alignmentFiles):
    '''takes a list of alignment file names.
       returns a concatenated alignment sorted by record id, with any taxa with missing data removed'''
    alignments = [AlignIO.read(file, get_alignment_type(file)) for file in alignmentFiles]
    
    numberOfTaxa = len(alignments[0])
    for alignment in alignments:
        if len(alignment) != numberOfTaxa:
            print("The alignments do not have equal numbers of taxa!")
            exit()
    
    for alignment in alignments: #make sure they're all sorted by name so it's all the same
        alignment.sort()
    
    
    #makes a set of all of the taxa missing in any alignment
    for alignment in alignments:
        alignment = replace_missing(alignment)
            
    concatenated = alignments[0]
    
    try:
        for alignment in alignments[1:]:
            concatenated += alignment
        
    except: # TODO find out what the actual exception is
        print('could not concatenate the alignments')
    
    
    return concatenated


def make_directories(dirList=[]):
    if dirList == []:
        dirList = ['accD', 'atp', 'ccsA', 'cemA', 'clpP', 'infA', 'matK', 'misc', 'ndh', 'pet', 'psa', 'psb', 'rbcL',  'rpl', 'rpo', 'rps', 'ycf3', 'ycf4']
    for name in dirList:
        if not os.path.isdir(name):
            os.mkdir(name)
    


def main():
    alignmentFiles = glob('*.fasta')
    geneType = 'plastid'
    if len(argv) > 1:
        geneType = argv[1]
    
    geneGroups = {}
    for file in alignmentFiles:
        groupNames = check_gene_group(get_gene_name(file), geneDict=geneType)
        for groupName in groupNames:
            if(groupName not in geneGroups.keys()):
                geneGroups[groupName] = [file]
            else:
                geneGroups[groupName].append(file)
    
    make_directories(geneGroups.keys())  # makes all of the directories needed
    
    for group, files in geneGroups.items():
        outFileName = group + '_' + '_'.join(files[0].split('_')[1:]) # makes the out file name based on the first file
        outFileName = os.path.join(group, outFileName)
        concatenatedAlignment = concat(files)
        AlignIO.write(concatenatedAlignment, handle=outFileName, format=get_alignment_type(alignmentFiles[0]))
    


if __name__ == '__main__':
    main()