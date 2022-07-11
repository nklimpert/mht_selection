# run this from base directory
# deletes all of the alignments in each taxa subdirectory of the alignments folder
# EXPECTS A LIST OF GENES IN THE ALIGNMENTS FOLDERS

dirnames=$1 # this will be the input

while read subdir; do #this will loop through each subdirectory
	echo $subdir
    origin="./alignments/$subdir" #should be the subdidrectory name



    files="$origin/*.fasta" #works
    regex="(\w+)_[a-zA-Z]+_.+" #works

    #HAS TO BE IN THE ALIGNMENTS FOLDER
    hits="$origin/genelist_$subdir.txt" # Give it the genelist for what genes to look for
    

    for f in $files
    do
        if [[ $f =~ $regex ]]
        then
            name="${BASH_REMATCH[1]}"

        
            if grep -Fxq "${name}" "$hits"
            then
                echo "$f should stay"
            else
                #echo "$f should be deleted"
                rm $f
            fi
            
        
        fi
    done 

done < $dirnames