#!/bin/bash
#
g++ -c -O3 md.cpp >& compiler.txt
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
mv a.out md_O3
./md_O3 < md_input.txt > md_O3_output.txt
rm md_O3
#
echo "Output written to md_O3_output.txt"
