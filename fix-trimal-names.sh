#!/bin/bash

parallel "sed -i'' 's/ [0-9]* bp//g' {}" ::: *.fasta