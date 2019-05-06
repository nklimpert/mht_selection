#!/bin/bash

parallel "cd {//} && python ../../scripts/filter.py {/} taxa_*.txt" ::: alignments/*/*.fasta
