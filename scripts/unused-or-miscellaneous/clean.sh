#/bin/zsh

directory=$1

if [ $directory == "alignments" ]; then
	for dir in $directory/*/*; do
		find $dir -type f -name '*.fasta' # -exec rm -rf {} \;
		find $dir -type d -name "*" # -exec rm -rf {} \;
		find $dir -type f -name "*.txt" # -exec rm -rf {} \;
	done
elif [ $directory == "paml" ]; then
	for dir in $directory/*/*; do
		find $dir -type f -name '*.fasta' # -exec rm -rf {} \;

	done
elif [ $directory == "relax" ]; then
	for dir in $directory/*/*; do
		find $dir -type f -name '*.fasta' # -exec rm -rf {} \;

	done
fi