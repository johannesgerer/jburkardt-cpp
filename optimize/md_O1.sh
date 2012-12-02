#!/bin/bash
#
g++ -c -O1 md.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md.cpp"
  exit
fi
rm compiler.txt
#
g++ md.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o"
  exit
fi
rm md.o
#
mv a.out md_O1
./md_O1 < md_input.txt > md_O1_output.txt
rm md_O1
#
echo "Output written to md_O1_output.txt"
