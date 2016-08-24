#!/bin/bash
for i in $(ls -1 *.afa); do
        g=$(wc -l $i | awk '{print $1}')
        if [ "$g" -lt 7 ]; then
                mv $i ./lessthanfour/$i
        fi
done