#run this from base mycoheterotroph-selection directory
#copies all of the genes from the alignments subfolders to the relax subfolders


dirnames=$1 # this will be the input

while read subdir; do #this will loop through each subdirectory

    origin="./alignments/$subdir" #should be the subdidrectory name
    target="./relax/$subdir"



    files="$origin/*.fasta" #works
    regex="(\w+)_[a-zA-Z]+_.+" #works
    hits="$target/genelist_$subdir.txt" # Give it the genelist for what genes to look for works


    for f in $files
    do
        if [[ $f =~ $regex ]]
        then
            name="${BASH_REMATCH[1]}"
            if grep -Fxq "${name}" "$hits"
            then

                echo "$f"
                genedir="$target/${name}/"
                #echo $genedir
                #TODO: uncomment the next line to make it copy file to directory
                cp $f "$genedir"
            
            fi
        
        else
            echo "$f doesn't match" >&2 # this could get noisy if there are a lot of non-matching files
        fi
    done 

done < $dirnames