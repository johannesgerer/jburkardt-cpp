#!/bin/bash
#
g++ -c tensor_example1.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tensor_example1.cpp"
  exit
fi
rm compiler.txt
#
g++ tensor_example1.o
if [ $? -ne 0 ]; then
  echo "Errors linking tensor_example1.o."
  exit
fi
#
rm tensor_example1.o
#
mv a.out tensor_example1
./tensor_example1 > tensor_example1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tensor_example1."
  exit
fi
rm tensor_example1
#
echo "Program output written to tensor_example1_output.txt"
