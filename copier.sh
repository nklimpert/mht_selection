while read lineage;
do
    cp ./petrosaviids/*.fasta ./$lineage/
    
done < ../lineages.txt