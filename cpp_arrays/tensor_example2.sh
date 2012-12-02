#!/bin/bash
#
g++ -c tensor_example2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tensor_example2.cpp"
  exit
fi
rm compiler.txt
#
g++ tensor_example2.o
if [ $? -ne 0 ]; then
  echo "Errors linking tensor_example2.o."
  exit
fi
#
rm tensor_example2.o
#
mv a.out tensor_example2
./tensor_example2 > tensor_example2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tensor_example2."
  exit
fi
rm tensor_example2
#
echo "Program output written to tensor_example2_output.txt"
