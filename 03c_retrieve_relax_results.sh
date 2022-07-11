#!/bin/bash

cd relax/
for clade in */; do
	cd $clade
	for genedir in */; do
		cd $genedir
		echo 'Retrieving results from' $genedir 'in' $clade
		echo *.json
		python ../../../scripts/retrieve_relax_results.py *.RELAX.json
		cd ..
	done
	cd ..
done

echo 'clade,gene,p-value,K-value' > combined_results.csv
find */*/*.csv -exec cat {} \; >> combined_results.csv