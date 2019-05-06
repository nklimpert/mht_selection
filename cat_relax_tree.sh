#

for dir in relax/*/; do
    cd $dir
    for subdir in */; do
        cd $subdir
        cat *.tre >> *.fasta
        cd ..
    done
    cd ../../
done