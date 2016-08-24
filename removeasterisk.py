#!/usr/bin/python
import re
from pdb import set_trace as bp
from Bio import SeqIO

# for i in range(1, 17261):
filenames = [line.strip() for line in open('filenames.txt', 'r')]
# bp()
for i in range(len(filenames)):
	sequences = list(SeqIO.parse('./' + filenames[i], "fasta"))
	for f in range(len(sequences)):
		if sequences[f].seq[-1] == '*':
			sequences[f].seq = sequences[f].seq[:-1]
			print('removing * from ' + str(sequences[f].id))
	print('writing ' + str(filenames[i]))
	SeqIO.write(sequences, './noast/' + filenames[i], 'fasta')


