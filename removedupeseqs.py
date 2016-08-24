#!/usr/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import re
import subprocess
from Bio.SeqIO import FastaIO
from pdb import set_trace as bp
import sys

# filenames = [line.strip() for line in open("./filenames.txt", 'r')] #grab taxa file names #list of files
# filenames = ["newmcldump/i14seqs/filenames.txt"] #grab taxa file names #single file

filenames = str(sys.argv[1])
log = open('./logs/' + filenames + '.log', 'a')

filecount = 0
# for loo in range(len(filenames)):
records = list(SeqIO.parse('./' + filenames, 'fasta'))
nodupes = 0
recstodel = []
matchindexno = []
rlength = range(len(records))
count = 0
rlength = range(len(records))
while count < len(records):
	overlap = False
	endcount = 0
	# matchindexno = []
	ex = records[count].seq
	exid = records[count].id
	for bla in range(count, len(records)):
		taxa = re.split("/", exid)
		# bp()
		if bla + 1 < len(rlength) and ex == records[bla + 1].seq:
			taxa1 = re.split("/", records[bla + 1].id)
			if bla + 1 not in matchindexno and taxa[0] == taxa1[0]:
				matchindexno.append(bla + 1)
				log.write("\nFOUND MATCH: " + str(records[bla + 1].id) + "\n")
			# matchindexno.append(bla + 1)
			# print matchindexno
			# bp()
			endcount += 1
			# break
	count += 1
log.write(str(matchindexno))
print("\n\n")
recstodel = matchindexno
log.write("Deleting " + str(len(recstodel)) + " records")
recstodel.sort()
print recstodel
# print recstodel
if recstodel != []:
	for x in range(len(recstodel)):
		# print("deleting record " + str(recstodel[x]))
		# print(recstodel)
		del records[recstodel[x]]
		recstodel = [f - 1 for f in recstodel]
else:
	pass
# print(len(records))
SeqIO.write(records, './cleaned/' + filenames, 'fasta')
log.close()


