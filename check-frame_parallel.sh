#!/bin/bash

parallel "cd {//} && sed 's/-/N/g' {/} > temp_{/} && python ../../scripts/check-frame.py temp_{/} >> issues.txt && rm temp_{/}" ::: alignments/*/*.fasta