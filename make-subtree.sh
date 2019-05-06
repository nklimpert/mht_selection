#!/bin/bash

lineage_list=$1

while read lineage; do
	for taxa_list in taxa_$lineage.txt; do
		python ../scripts/make-subtree.py petrosaviids.tre $taxa_list $lineage
	done
done < $lineage_list