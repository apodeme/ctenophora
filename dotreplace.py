#!/usr/bin/python
#replaces dots in gene names (for ALE)
import sys
filenames = str(sys.argv[1])
filedata = open(filenames, 'r')
taxa = list(filedata.read())

for i in range(len(taxa)):
	# print taxa[i]
	if taxa[i] == '.':
		if taxa[i - 2] == ':' or taxa[i - 3] == ':' or taxa[i - 4] == ':' or taxa[i - 5] == ':':
			pass
		else:
			taxa[i] = 'DOTHERE'
strtaxa = "".join(taxa)

newfile = open('./nodots2/' + filenames,'w')
newfile.write(str(strtaxa))
filedata.close
newfile.close

# parallel --joblog paradots2.log './nodots.py {}' ::: *.afa.ufboot