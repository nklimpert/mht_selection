#!usr/bin/env python

from scipy.stats import chi2
import os
from sys import argv


class TestResults:
    """results from individual tests of paml"""
    def __init__(self, clade='', gene='', model='', lnL=0, *args):
        self.clade = clade
        self.gene = gene
        self.model = model
        self.model_name = model.split('.')[0]
        self.lnL = float(lnL)

        tests = {'M0': 'branch',
                 'b_free': 'branch',
                 'bsA1': 'bsA',
                 'bsA': 'bsA',
                 'M3': 'cmD',
                 'bsD': 'cmD',
                 'XX': 'cmC',
                 'bsC': 'cmC'
                 }
        self.test = tests[self.model_name]

        extra_parameters = [float(a) for a in args]
        model_types = {'M3': lambda l:
                       {'Proportions': l[0::2], 'Omegas': l[1::2]},
                       'M0': lambda l:
                       {'Omegas': l},
                       'bsD': lambda l:
                       {'Proportions': l[0::3], 'Background_omegas': l[1::3],
                        'Foreground_omegas': l[2::3]},
                       'b_free': lambda l:
                       {'Foreground_omegas': l[0::2], 'Background_omegas': l[1::2]},
                       'XX': lambda l:
                       {'Proportions': l[0::2], 'Omegas': l[1::2]},
                       'bsC': lambda l:
                       {'Proportions': l[0::3], 'Background_omegas': l[1::3],
                        'Foreground_omegas': l[2::3]}
                       }
        # this is where the omegas are stored
        self.parameters = model_types[self.model_name](extra_parameters)

    def __str__(self):
        return str(self.clade + '_' + self.gene + '_' + self.model)

    def __repr__(self):
        return str(self.clade + '_' + self.gene + '_' + self.model)

    def __lt__(self, other):
        return self.lnL < other.lnL

    def __le__(self, other):
        return self.lnL <= other.lnL

    def __eq__(self, other):
        return self.lnL == other.lnL

    def __ne__(self, other):
        return self.lnL != other.lnL

    def __gt__(self, other):
        return self.lnL > other.lnL

    def __ge__(self, other):
        return self.lnL >= other.lnL


class PamlResults:
    '''encapsulate all the different paml tests run for a clade+gene combo'''

    def likelihood_ratio(self, llmin, llmax):
        ''' returns the likelihood ratio of the two log likelihoods'''
        return (2*abs(llmax-llmin))

    def chi_dist(self, LRT):
        ''' returns the p-value based on the chi squared distribution on the likelihood ratio test
            with one degree of freedom'''
        return chi2.sf(LRT, 1)

    def __init__(self, clade, gene, *args):
        self.clade = clade
        self.gene = gene

        # args is going to be a series of lists, each with two elements
        # these are dependent upon what kinds of models are being tested - each kind of test will
        # have two models
        # one model for the null and one for the test

        self.tests = {}  # key will be the kind of test it was, item will be [H0, H1] models
        for modelPair in args:
            testName = modelPair[0].test
            self.tests[testName] = modelPair

        # next, make a dictionary of the likelihood ratio tests

        self.LRTs = {}
        for testName, modelPair in self.tests.items():
            self.LRTs[testName] = self.likelihood_ratio(modelPair[0].lnL, modelPair[1].lnL)

        # next, make a dictionary of the p-values

        self.p_values = {}
        for testName, likelihood in self.LRTs.items():
            self.p_values[testName] = self.chi_dist(likelihood)

        # initialize all of the fdr values to -1
        self.fdr_values = {}
        for testName in self.tests.keys():
            self.fdr_values[testName] = -1

        self.significance = {}
        for testName in self.tests.keys():
            self.significance[testName] = None

    def __str__(self):
        values = [self.clade, self.gene]
        for testName in self.tests.keys():
            values.extend([self.LRTs[testName], self.p_values[testName], self.fdr_values[testName],
                           self.significance[testName]])
        return ','.join(map(lambda x: str(x), values))

    def __repr__(self):
        values = [self.clade, self.gene]
        for testName in self.tests.keys():
            values.extend([self.LRTs[testName], self.p_values[testName], self.fdr_values[testName],
                           self.significance[testName]])
        return ','.join(map(lambda x: str(x), values))


def fdr_correction(results):
    ''' does false discovery rate correction with Benjamin-Hochberg procedure on a list of PamlResults
        returns a dictionary where the name is the key and the value is a boolean
        (True if significant)'''
    for testName in results[0].tests.keys():
        sorted_results = sorted(results, key=lambda x: x.p_values[testName])
        for i in range(len(sorted_results)):
            fdr_p = 0.05*(float(i+1)/len(sorted_results))
            sorted_results[i].fdr_values[testName] = fdr_p
            sorted_results[i].significance[testName] = sorted_results[i].p_values[testName] \
                < sorted_results[i].fdr_values[testName]


def file_reader(resultsFile):
    all_results = []
    with open(resultsFile, 'r') as csvfile:
        for line in csvfile:
            # get rid of all of the empty cells
            row = [x for x in filter(lambda x: x != '', line.rstrip().split(','))]
            if 'gene' not in row:  # make it skip the first row only if it's a title row
                all_results.append(TestResults(*row))
    return all_results


def main():

    resultsFile = argv[1]

    all_results = file_reader(resultsFile)  # list of results

    cladeDict = {}
    for result in all_results:
        if result.clade not in cladeDict.keys():
            cladeDict[result.clade] = {result.gene: {}}
        if result.gene not in cladeDict[result.clade].keys():
            cladeDict[result.clade][result.gene] = {}
        if result.test not in cladeDict[result.clade][result.gene].keys():
            cladeDict[result.clade][result.gene][result.test] = []
        cladeDict[result.clade][result.gene][result.test].append(result)

    all_PamlResults = {}

    for clade in cladeDict.keys():
        all_PamlResults[clade] = []
        for gene in cladeDict[clade].keys():
            #TODO make it so that it puts all of the tests together
            listOfTests = [v for v in cladeDict[clade][gene].values()]

            pResult = PamlResults(clade, gene, *listOfTests)
            all_PamlResults[clade].append(pResult)
        fdr_correction(all_PamlResults[clade])

    with open(os.path.splitext(resultsFile)[0]+'_corrected.csv', 'w+') as outFile:
        header = "clade,gene"
        # look at the first result
        for testName in list(all_PamlResults.values())[0][0].tests.keys():
            header += ",{0} LRT, {0} p-value, {0} FDR, {0} significant?".format(testName)
        header += "\n"
        outFile.write(header)
        for genes in all_PamlResults.values():
            for gene in genes:
                outFile.write(str(gene)+'\n')


if __name__ == '__main__':
    main()
