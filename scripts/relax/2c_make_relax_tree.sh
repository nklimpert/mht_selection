#!/bin/bash

cd relax/
for clade in */; do
	cd $clade
	for gene in */; do
		cd $gene
		echo "Creating a tree for" "${gene:0:-1}" in $clade
		python ../../../scripts/make_relax_tree.py ../petrosaviids_"${clade:0:-1}".tre "${gene:0:-1}"_"${clade:0:-1}"_2018-*.fasta ../test_taxa_"${clade:0:-1}".txt
		cd ../
	done
	cd ../
done