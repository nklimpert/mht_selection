#!/bin/bash

source activate ete3
parallel -j3 -u --joblog paml.log "cd {} && python ../../../scripts/paml.py *.fasta ../*.tre branch ../test_taxa_*.txt" ::: paml/*/*/
