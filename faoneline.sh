#!/bin/bash
#convert to one line fasta

#convert all files in txt from multiline fasta to single line
count=1
numfiles=$(($(wc -l ./taxafilenames.txt | awk '{print $1}') + 1)) 
echo $numfiles
while [ $count -le $numfiles ]; do
	curfile=$(sed "${count}q;d" ./taxafilenames.txt)
	curfile=$(echo "$curfile")
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ./$curfile > ./oneline/$curfile
	count=$(($count + 1))
done



#remove first blank line produced by awk with sed
count=1
numfiles=$(($(wc -l ./taxafilenames.txt | awk '{print $1}') + 1)) 
echo $numfiles
while [ $count -le $numfiles ]; do
	curfile=$(sed "${count}q;d" ./taxafilenames.txt)
	curfile2=$(sed "${count}q;d" ./taxafilenames.txt)
	curfile=$(echo "$curfile")
	sed 1d ./oneline/$curfile > ./oneline/wolineone/$curfile2
	count=$(($count + 1))
done

