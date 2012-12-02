#!/bin/bash
#
g++ -c -g md.cpp >& compiler.txt
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
mv a.out md
./md < md_input.txt > md_output.txt
rm md
#
echo "Output written to md_output.txt"
