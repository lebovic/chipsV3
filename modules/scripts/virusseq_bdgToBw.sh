#!/bin/bash
#CHECKS to see if the input file (arg 1) is empty; IF it is, then touches the
#output file (arg 3),
#otherwise runs bedGraphToBigWig {input} {assembly len- arg2} {output}

#check if input empty
if [ ! -s $1 ]; then
    #echo "File is empty"
    touch $3
else
    #Run bedGraphToBigWig
    bedGraphToBigWig $1 $2 $3
fi
