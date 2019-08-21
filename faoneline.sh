#!/bin/bash
#convert to one line fasta
#convert all files in txt from multiline fasta to single line
mkdir ./oneline

count=1
numfiles=$(($(ls -l ./*.fa | wc -l | awk '{print $1}') - 1)) 
echo $numfiles
while [ $count -le $numfiles ]; do
	curfile=$(ls -l | awk '{print $9}' | sed "${count}q;d")
	curfile=$(echo "$curfile")
	echo $curfile
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ./$curfile > ./oneline/$curfile
	count=$(($count + 1))
done


mkdir ./oneline/wolineone
#remove first blank line produced by awk with sed
count=1
numfiles=$(($(ls -l ./*.fa | wc -l | awk '{print $1}') - 1))  
echo $numfiles
while [ $count -le $numfiles ]; do
	curfile=$(ls -l | awk '{print $9}' | sed "${count}q;d")
	curfile2=$(ls -l | awk '{print $9}' | sed "${count}q;d")	
	curfile=$(echo "$curfile")
	sed 1d ./oneline/$curfile > ./oneline/wolineone/$curfile2
	count=$(($count + 1))
done

mkdir multi.fa
mv ./*.fa multi.fa

mv ./oneline/wolineone/*.fa ./
rm -r ./oneline
