'''take two lists, create a list that contains all the elements unique to the first list
    Usage: difference_list.py <fullList> <file2>
'''

from sys import argv

file1 = argv[1]

file2 = argv[2]

with open(file1, 'r') as inFile:
    file1_set = {line.strip() for line in inFile.readlines()}

with open(file2, 'r') as inFile:
    file2_set = {line.strip() for line in inFile.readlines()}

diff = file1_set - file2_set

print('\n'.join(diff))
