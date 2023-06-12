
for dir in relax/*/; do
    cd $dir
    for subdir in */; do
        cd $subdir
        echo "Making a bf for" *.fasta "in dir" $subdir
        python ../../../scripts/make_bf.py *.fasta
        cd ..
    done
    cd ../../
done