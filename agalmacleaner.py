#!/usr/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import re
import subprocess
from Bio.SeqIO import FastaIO
# from pdb import set_trace as bp

# filenames = [line.strip() for line in open("taxafilenames.txt", 'r')] #grab taxa file names #list of files
filenames = ["testseqs.fa"] #grab taxa file names #single file

def all_same(items):
    return all(x == items[0] for x in items)

filecount = 0
for loo in range(len(filenames)):
	records = list(SeqIO.parse('./alltaxaaas/' + filenames[loo], 'fasta'))
	nodupes = 0
	recstodel = []
	matchindexno = []
	rlength = range(len(records))
	count = 0
	rlength = range(len(records))
	while count < len(records):
		overlap = False
		endcount = 0
		matchindexno = []
		ex = records[count].id
		for bla in range(count, len(records)):
			if bla + 1 < len(rlength) and ex == records[bla + 1].id:
				matchindexno.append(bla)
				endcount += 1
				print("\nFOUND MATCH: " + str(records[bla].id))
			elif endcount != 0:
				print("FOUND MATCH: " + str(records[bla].id))
				matchindexno.append(bla)
				matchindexnoadd1 = [x+1 for x in matchindexno]
				print("INDEX POSITIONS: " + str(matchindexnoadd1))
				endcount = 0 #endcount must be set back to zero here!
				# print(matchindexno)
				# for i in matchindexno:
				# 	print records[i].id
			else:
				endcount = 2
				break

		if matchindexno != []:
			signsl = []
			for par in range(len(matchindexno)): #parse signs and add to list
				split = re.split(' ', records[matchindexno[par]].description)
				porm = re.search('\(\D\)', split[1])
				signsl.append(porm.group(0)[1])
			print signsl

			if all_same(signsl) == True:
				print("SIGNS ALL SAME")
				framenos = []
				for fra in range(len(matchindexno)): #parse frame ranges
					split = re.split(' ', records[matchindexno[fra]].description)
					frame = re.search('\d*-\d*', split[1])
					framerange = re.split('-', frame.group(0))
					framerange = [int(x) for x in framerange]
					framerange.append("Pos.:")
					framerange.append(matchindexno[fra]) #add list index number to ranges for after sorting
					framenos.append(framerange)
				print framenos
				framesorted = sorted(framenos, key=lambda x: x[0]) #sort ranges by start numbers
				w_it = 0
				overlap = False
				print('FRAMESORTED: ' + str(framesorted))
				print('len(framesorted) = ' + str(len(framesorted)))

				#CHECK IF SEQUENCES OVERLAP
				while w_it < len(framesorted) - 1 and overlap == False:
					print('w_it: ' + str(w_it))
					if framesorted[w_it][1] < framesorted[w_it + 1][0]:
						print('framesorted[w_it][1]: ' + str(framesorted[w_it][1]))
						print('framesorted[w_it + 1][0]: ' + str(framesorted[w_it + 1][0]))
						print('PASSING')
						# print('w_it: ' + str(w_it))
						pass
						# if w_it == len(framesorted) - 1
					else:
						overlap = True
						print('SETTING OVERLAP TRUE')
					w_it += 1
					# print('w_it: ' + str(w_it))

				if overlap == True: #if the sequences overlap in range, find biggest and throw away the rest
					print("SEQUENCES OVERLAP")
					totlen = []
					listtots = []
					print('SORTED FRAME LENGTHS: ' + str(framesorted))
					for i3 in range(len(framesorted)):
						totlen.append(int(framesorted[i3][1]) - int(framesorted[i3][0]))
						totlen.append('Pos.:')
						totlen.append(framesorted[i3][3])
						listtots.append(totlen)
						totlen = []
						print('LISTTOTS: ' + str(listtots))
					biggest = [0]
					for i4 in range(len(listtots)):
						if listtots[i4][0] > biggest[0]:
							biggest = [listtots[i4][0]]
							biggest.append(listtots[i4][1])
							biggest.append(listtots[i4][2])
					print("Biggest:" + str(biggest[0]))
					for yp in range(len(listtots)): #add tags for deletion
						if listtots[yp][2] == biggest[2]:
							listtots[yp].append("keep")
						else:
							listtots[yp].append("delete")
					print listtots
					delcount = 0
					for yf in range(len(listtots)): #delete smaller sequences
						if listtots[yf][3] == "delete":
							# print records[listtots[yf][2]]
							print("DELETING RECORD " + str(listtots[yf][2]) + ", ID: " + str(records[listtots[yf][2]].id))
							recstodel.append(listtots[yf][2])
							# records[listtots[yf][2]] = 0
							delcount += 1
							# print records[listtots[yf][2]]
					count += delcount + 1

				if overlap == False:
					if signsl[0] == '+': #if signs are positive, stick together in order specified by frame
						print('NO OVERLAP, SIGNS +, CONCATENATE')
						for gx in range(0, len(framesorted) - 1): #remove asterisks from all but the final sequence to be concatenated
							if records[framesorted[gx][3]].seq[-1] == '*':
								records[framesorted[gx][3]].seq = records[framesorted[gx][3]].seq[0:-1]
						for lx in range(1, len(framesorted)): #concatenate sequences
							records[framesorted[0][3]].seq = records[framesorted[0][3]].seq + records[framesorted[lx][3]].seq
						delcount = 0
						for px in range(1, len(framesorted)): #delete extra sequences
							print("DELETING RECORD " + str(framesorted[px][3] + 1) + ", ID: " + str(records[framesorted[px][3]].id))
							recstodel.append(framesorted[px][3])
							# records[framesorted[px][3]] = 0
							delcount += 1
						count += delcount + 1
					else: #if signs are negative, stick together in reverse order
						framesorted = sorted(framenos, key=lambda x: x[1], reverse=True) #sort ranges by end numbers, descending order
						print('NO OVERLAP, SIGNS -, CONCATENATE REVERSE')
						for gx in range(0, len(framesorted) - 1): #remove asterisks from all but the final sequence to be concatenated
							if records[framesorted[gx][3]].seq[-1] == '*':
								records[framesorted[gx][3]].seq = records[framesorted[gx][3]].seq[0:-1]
						for lx in range(1, len(framesorted)): #concatenate sequences
							records[framesorted[0][3]].seq = records[framesorted[0][3]].seq + records[framesorted[lx][3]].seq
						delcount = 0
						for px in range(1, len(framesorted)): #delete extra sequences
							print("DELETING RECORD " + str(framesorted[px][3] + 1) + ", ID: " + str(records[framesorted[px][3]].id))
							recstodel.append(framesorted[px][3])
							# records[framesorted[px][3]] = 0
							delcount += 1
						count += delcount + 1

			if all_same(signsl) == False:
				print("SIGNS NOT THE SAME")
				framenos = [] #repeat of previous code for selecting largest
				for fra in range(len(matchindexno)): #parse frame ranges
					split = re.split(' ', records[matchindexno[fra]].description)
					frame = re.search('\d*-\d*', split[1])
					framerange = re.split('-', frame.group(0))
					framerange = [int(x) for x in framerange]
					framerange.append("Pos.:")
					framerange.append(matchindexno[fra]) #add list index number to ranges for after sorting
					framenos.append(framerange)
				framesorted = sorted(framenos, key=lambda tup: tup[0]) #sort ranges by start numbers
					#find biggest range
				totlen = []
				listtots = []
				print('SORTED FRAME LENGTHS: ' + str(framesorted))
				for i3 in range(len(framesorted)):
					totlen.append(int(framesorted[i3][1]) - int(framesorted[i3][0]))
					totlen.append('Pos.:')
					totlen.append(framesorted[i3][3])
					listtots.append(totlen)
					totlen = []
					# print('LISTTOTS: ' + str(listtots))
				biggest = [0]
				for i4 in range(len(listtots)):
					if listtots[i4][0] > biggest[0]:
						biggest = [listtots[i4][0]]
						biggest.append(listtots[i4][1])
						biggest.append(listtots[i4][2])
				print("Biggest:" + str(biggest[0]))
				for yp in range(len(listtots)): #add tags for deletion
					if listtots[yp][2] == biggest[2]:
						listtots[yp].append("keep")
					else:
						listtots[yp].append("delete")
				print listtots
				delcount = 0
				for yf in range(len(listtots)): #delete smaller sequences
					if listtots[yf][3] == "delete":
						# print records[listtots[yf][2]] #debug: checking correct deletion
						print("DELETING RECORD " + str(listtots[yf][2] + 1) + ", ID: " + str(records[listtots[yf][2]].id))
						recstodel.append(listtots[yf][2])
						# records[listtots[yf][2]] = 0
						delcount += 1
						# print records[listtots[yf][2]] #debug: checking correct deletion
				count += delcount + 1
		else:
			count += 1

	print("\n\n")
	print("Deleting " + str(len(recstodel)) + " records")
	recstodel.sort()
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
	# SeqIO.write(records, './alltaxaaas/mod/' + filenames[loo] + '.mod', 'fasta')


	# print(str(nodupes) + "duplicate sequences detected")
