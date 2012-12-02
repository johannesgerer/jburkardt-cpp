#!/bin/bash
#
g++ -c sizes.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sizes.cpp"
  exit
fi
rm compiler.txt
#
g++ sizes.o
if [ $? -ne 0 ]; then
  echo "Errors loading sizes.o"
  exit
fi
rm sizes.o
#
mv a.out sizes
./sizes > sizes_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sizes"
  exit
fi
rm sizes
#
echo "Program output written to sizes_output.txt"
