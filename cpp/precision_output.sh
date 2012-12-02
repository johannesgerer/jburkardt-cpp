#!/bin/bash
#
g++ -c precision_output.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling precision_output.cpp"
  exit
fi
rm compiler.txt
#
g++ precision_output.o
if [ $? -ne 0 ]; then
  echo "Errors linking precision_output.o."
  exit
fi
#
rm precision_output.o
#
mv a.out precision_output
./precision_output > precision_output_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running precision_output."
  exit
fi
rm precision_output
#
echo "Program output written to precision_output_output.txt"
