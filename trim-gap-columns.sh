for dir in alignments/*/;
do
    cd $dir
    for alignment in *.fasta;
    do
    	echo "Working on" $alignment "in" $dir
        trimal -in $alignment -out $alignment.noallgaps -noallgaps
        sed 's/ [0-9]* bp//g' $alignment.noallgaps > $alignment
        rm $alignment.noallgaps
    done
    cd ../../
done