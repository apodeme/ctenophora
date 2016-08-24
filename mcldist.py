#!/usr/bin/python
import re
# from pdb import set_trace as bp
from Bio import SeqIO
import sys
filenames = str(sys.argv[1])
# filenames =

#load in MCL clusters

mclfile_handle = open("./" + str(filenames))

# check number of fields in each row
count = 1
output_filename = str(filenames) + "_fdist.txt"
output_handle = open(output_filename, "a")
output_handle.write("row fieldno\n")
output_handle.close()
for row in mclfile_handle:
	seqnames = re.split("\t", row.rstrip())
	output_handle = open(output_filename, "a")
	output_handle.write(str(count) + " " + str(len(seqnames)) + "\n")
	output_handle.close()
	count += 1

